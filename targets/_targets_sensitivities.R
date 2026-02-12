suppressPackageStartupMessages(suppressWarnings({
  library(targets)
  library(crew)
  library(magrittr)
  library(qs)
  library(tidyverse)
}))

tar_option_set(
  format = "qs", 
  controller = crew::crew_controller_local(
    workers = 8, 
    seconds_idle = 120
  ),
  workspace_on_error = TRUE,
  garbage_collection = 10
)

tar_source(
  files = c(
    here::here("src/dirs.R"),
    here::here("src/functions.R"),
    here::here("src/model_functions.R")
  )
)

# Dirs
output_path <- here() %>% file.path("outputs")
output_species_data_path <- file.path(output_path, "species_data")

# Globals
farm_chunk_size <- 10
inds_per_farm <- 2500

list(
  # Load previously saved data --------------------------------------------------------------------------------
  tar_target(species_params_file, file.path(output_species_data_path, "species_params.qs"), format = "file"),
  tar_target(pop_params_file, file.path(output_species_data_path, "pop_params.qs"), format = "file"),
  tar_target(feed_params_file, file.path(output_species_data_path, "feed_params_PD.qs"), format = "file"),
  tar_target(farm_temporal_data_file, file.path(input_farm_coords_path, "farm_temporal_data.qs"), format = "file"),
  tar_target(farm_static_data_file, file.path(input_farm_coords_path, "farm_static_data_w_stocking.qs"), format = "file"),

  tar_target(
    farm_IDs,
    command = {
      sample_size <- farm_static_data_file %>% 
        qread() %>% pull(farm_id) %>% unique() %>% length()
      sample_size <- round(sample_size * 0.1, 0) # 10% of total

      farm_static_data_file %>% 
        qread() %>% 
        distinct(farm_id) %>% 
        slice_sample(n = sample_size) %>% 
        pull(farm_id)
    }
  ),

  tar_target(
    farm_IDs_chunked,
    split(farm_IDs, ceiling(seq_along(farm_IDs)/farm_chunk_size))
  ),

  tar_target(feed_params_ref, qread(feed_params_file)[["plant_dominant_biomar"]]),

  tar_target(
    test_reference_feed,
    command = {
      fi <- sample(farm_IDs, 1)

      st_data <- farm_static_data %>% 
        filter(farm_ID == fi)
      ts_data <- farm_temporal_data_file %>% 
        qread() %>% 
        filter(farm_id == fi) %>% 
        arrange(day) %>% 
        filter(day >= st_data$t_start & day <= st_data$t_end)

      fg <- fish_growth(
        pop_params = all_params,
        species_params = all_params,
        water_temp = ts_data$temp_c,
        feed_params = feed_params_ref,
        times = c(t_start = st_data$t_start, t_end = st_data$t_end, dt = 1),
        init_weight = all_params["meanW"],
        ingmax = all_params["meanImax"]
      )
      fg %>%
        as.data.frame() %>%
        filter(!is.na(weight)) %>% 
        slice_tail(n = 1) %>%
        select(weight) %>%
        mutate(
          farm_ID = fi,
          weight = weight %>% 
            units::set_units("g") %>% 
            units::drop_units()
        )
    }
  ),

# Prepare parameters and farm data --------------------------------------------------------------------------------
  tar_target(
    farm_static_data,
    command = {
      farm_static_data_file %>% 
        qread() %>% 
        sf::st_drop_geometry() %>% 
        rename(farm_ID = farm_id) %>% 
        filter(farm_ID %in% farm_IDs_chunked[[1]])
    },
    pattern = farm_IDs_chunked
  ),

  tar_target(
    farm_ts_data,
    command = {
      fsd <- farm_static_data %>% 
        filter(farm_ID %in% farm_IDs_chunked[[1]]) %>% 
        group_by(farm_ID) %>% 
        group_split()
      ftd <- farm_temporal_data_file %>% 
        qread() %>% 
        rename(farm_ID = farm_id) %>% 
        arrange(farm_ID, day) %>% 
        filter(farm_ID %in% farm_IDs_chunked[[1]]) %>% 
        group_by(farm_ID) %>% 
        group_split()

      # For testing
      # static <- fsd[[1]]
      # ts <- ftd[[1]]

      # Add population from stocking n
      purrr::map2_dfr(
        fsd, 
        ftd, 
        function(static, ts) {
          cbind(
            ts %>% 
              filter(day >= static$t_start & day <= static$t_end),
            data.frame(Npop = generate_pop(
                harvest_n = static$harvest_n,
                mort = all_params['mortmyt'],
                times = c(t_start = static$t_start, t_end = static$t_end, dt = 1)
            ))
          )
        }
      )
    },
    pattern = farm_IDs_chunked
  ),

  # Get names of params to be adjusted
  tar_target(all_params, c(qread(species_params_file), qread(pop_params_file))),
  tar_target(param_names_pop, c("meanW", "deltaW", "meanImax", "deltaImax", "overFmean", "overFdelta")),
  tar_target(param_names_spec, names(all_params)[!names(all_params) %in% param_names_pop]),

  tar_target(factors, c(0.9, 1, 1.1)),
  tar_target(
    sens_run_spec, 
    command = {
      # Split static and timeseries data by farm
      ts_data <- farm_ts_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      static_data <- farm_static_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      # Adjust one parameter
      adj_params <- all_params
      adj_params[param_names_spec] <- adj_params[param_names_spec] * factors

      # For testing
      # ts <- ts_data[[1]]
      # static <- static_data[[1]]

      # Run uniform farm growth for all farms
      purrr::map2_dfr(ts_data, static_data, function(ts, static) {
        fg <- uni_farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = feed_params_ref,
          times = c(t_start = static$t_start, t_end = static$t_end, dt = 1),
          N_pop = ts$Npop
        )
        fg %>% 
          purrr::map_dfr(function(mat) {
            mat %>% 
              as.data.frame() %>% 
              rename(t = V1, mean = V2) %>% 
              mutate(
                farm_ID = as.integer(static$farm_ID),
                t = as.integer(t)
              )
          }, .id = "measure") %>% 
          mutate(measure = as.factor(measure))
      }) %>% 
        # Label by param adjusted and adjustment factor
        mutate(
          adj_param = as.factor(param_names_spec),
          factor = factors
        ) %>% 
        group_by(farm_ID) %>% 
        mutate(slice_t = max(t)) %>% 
        ungroup() %>% 
        mutate(slice_t = case_when(measure == "weight_stat" ~ slice_t, T ~ slice_t-1)) %>% 
        filter(t == slice_t) %>%
        select(-t, -slice_t)
    },
    pattern = cross(param_names_spec, factors)
  ),

  tar_target(
    sens_results_spec,
    command = {
      sens <- sens_run_spec %>%
        pivot_wider(
          names_from = "factor", 
          values_from = "mean", 
          names_prefix = "fact_"
        ) %>% 
        mutate(sensitivity = (fact_1.1-fact_0.9)/(0.2*fact_1)) 
      
      sens %>% 
        group_by(adj_param, measure) %>%
        reframe(
          mean_sens = meanna(sensitivity),
          sd_sens = sdna(sensitivity)
        )
    }
  ),

  tar_target(
    sens_run_pop_lo, 
    command = {
      # Split static and timeseries data by farm
      ts_data <- farm_ts_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      static_data <- farm_static_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      # Adjust one parameter
      adj_params <- all_params
      adj_params[param_names_pop] <- adj_params[param_names_pop] * factors[1]

      # For testing
      # ts <- ts_data[[1]]
      # static <- static_data[[1]]
      # tmp <- 

      # Run growth for all farms
      purrr::map2_dfr(ts_data, static_data, function(ts, static) {
        fg <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = feed_params_ref,
          times = c(t_start = static$t_start, t_end = static$t_end, dt = 1),
          N_pop = ts$Npop,
          nruns = inds_per_farm
        )
        fg %>% 
          purrr::map_dfr(function(mat) {
            mat %>% 
              as.data.frame() %>% 
              select(-V3) %>% 
              rename(t = V1, mean = V2) %>% 
              mutate(
                farm_ID = as.integer(static$farm_ID),
                t = as.integer(t)
              )
          }, 
          .id = "measure"
        ) %>% 
          mutate(measure = as.factor(measure))
      } 
      # , .progress = T
      ) %>% 
        # Label by param adjusted and adjustment factor
        mutate(
          adj_param = as.factor(param_names_pop),
          factor = factors[1]
        ) %>% 
        group_by(farm_ID) %>% 
        mutate(slice_t = max(t)) %>% 
        ungroup() %>% 
        mutate(slice_t = case_when(measure == "weight_stat" ~ slice_t, T ~ slice_t-1)) %>% 
        filter(t == slice_t) %>%
        select(-t, -slice_t)
    },
    pattern = cross(param_names_pop, map(farm_ts_data, farm_static_data))
  ),

  tar_target(
    sens_run_pop_mi, 
    command = {
      # Split static and timeseries data by farm
      ts_data <- farm_ts_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      static_data <- farm_static_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      # Adjust one parameter
      adj_params <- all_params
      adj_params[param_names_pop] <- adj_params[param_names_pop] * factors[2]

      # For testing
      # ts <- ts_data[[1]]
      # static <- static_data[[1]]

      # Run growth for all farms
      purrr::map2_dfr(ts_data, static_data, function(ts, static) {
        fg <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = feed_params_ref,
          times = c(t_start = static$t_start, t_end = static$t_end, dt = 1),
          N_pop = ts$Npop,
          nruns = inds_per_farm
        )
        fg %>% 
          purrr::map_dfr(function(mat) {
            mat %>% 
              as.data.frame() %>% 
              select(-V3) %>% 
              rename(t = V1, mean = V2) %>% 
              mutate(
                farm_ID = as.integer(static$farm_ID),
                t = as.integer(t)
              )
          }, .id = "measure") %>% 
          mutate(measure = as.factor(measure))
      }) %>% 
        # Label by param adjusted and adjustment factor
        mutate(
          adj_param = as.factor(param_names_pop),
          factor = factors[2]
        ) %>% 
        group_by(farm_ID) %>% 
        mutate(slice_t = max(t)) %>% 
        ungroup() %>% 
        mutate(slice_t = case_when(measure == "weight_stat" ~ slice_t, T ~ slice_t-1)) %>% 
        filter(t == slice_t) %>%
        select(-t, -slice_t)
    },
    pattern = cross(param_names_pop, map(farm_ts_data, farm_static_data))
  ),

  tar_target(
    sens_run_pop_hi, 
    command = {
      # Split static and timeseries data by farm
      ts_data <- farm_ts_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      static_data <- farm_static_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      # Adjust one parameter
      adj_params <- all_params
      adj_params[param_names_pop] <- adj_params[param_names_pop] * factors[3]

      # For testing
      # ts <- ts_data[[1]]
      # static <- static_data[[1]]

      # Run growth for all farms
      purrr::map2_dfr(ts_data, static_data, function(ts, static) {
        fg <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = feed_params_ref,
          times = c(t_start = static$t_start, t_end = static$t_end, dt = 1),
          N_pop = ts$Npop,
          nruns = inds_per_farm
        )
        fg %>% 
          purrr::map_dfr(function(mat) {
            mat %>% 
              as.data.frame() %>% 
              select(-V3) %>% 
              rename(t = V1, mean = V2) %>% 
              mutate(
                farm_ID = as.integer(static$farm_ID),
                t = as.integer(t)
              )
          }, .id = "measure") %>% 
          mutate(measure = as.factor(measure))
      }) %>% 
        # Label by param adjusted and adjustment factor
        mutate(
          adj_param = as.factor(param_names_pop),
          factor = factors[3]
        ) %>% 
        group_by(farm_ID) %>% 
        mutate(slice_t = max(t)) %>% 
        ungroup() %>% 
        mutate(slice_t = case_when(measure == "weight_stat" ~ slice_t, T ~ slice_t-1)) %>% 
        filter(t == slice_t) %>%
        select(-t, -slice_t)
    },
    pattern = cross(param_names_pop, map(farm_ts_data, farm_static_data))
  ),

  tar_target(
    sens_results_pop,
    command = {
      sens <- rbind(
        sens_run_pop_lo,
        sens_run_pop_mi,
        sens_run_pop_hi
      )
      
      sens <- sens %>%
        pivot_wider(
          names_from = "factor", 
          values_from = "mean", 
          names_prefix = "fact_"
        ) %>% 
        mutate(sensitivity = (fact_1.1-fact_0.9)/(0.2*fact_1)) 
      
      sens %>% 
        group_by(adj_param, measure) %>%
        reframe(
          mean_sens = meanna(sensitivity),
          sd_sens = sdna(sensitivity)
        )
    }
  )
)
