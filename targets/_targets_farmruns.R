library(targets)
# library(crew)
library(magrittr)

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "tibble", "qs"), 
  format = "qs", 
  controller = crew::crew_controller_local(workers = 12, seconds_idle = 300),
  workspace_on_error = TRUE
)

tar_source(
  files = list.files("src", pattern = "\\.R$", full.names = TRUE) %>% 
    setdiff("src/map_templates.R")
)

list(
# Get previously saved data ---------------------------------------------------------------------------------------
  tar_target(farm_coords_file, file.path(output_farm_data_path, "farm_coords.qs"), format = "file"),
  tar_target(farm_ts_data_file, file.path(output_farm_data_path, "farm_ts_data.qs"), format = "file"),
  tar_target(species_params_file, file.path(output_species_data_path, "species_params.qs"), format = "file"),
  tar_target(pop_params_file, file.path(output_species_data_path, "pop_params.qs"), format = "file"),
  tar_target(feed_params_file, file.path(output_species_data_path, "feed_params.qs"), format = "file"),
  tar_target(farm_production_file, file.path(input_farm_coords_path, "atlantic_salmon_locations_w_temps.qs"), format = "file"),

  tar_target(farm_coords, qs::qread(farm_coords_file)),
  tar_target(
    farm_IDs,
    command = {
      qs::qread(farm_ts_data_file) %>% 
        filter(farm_ID != 2703) %>% # The only Russian farm
        distinct(farm_ID) %>% 
        pull(farm_ID) #%>% sample(112)
    }
  ),

  tar_target(
    farm_IDs_chunked,
    split(farm_IDs, ceiling(seq_along(farm_IDs)/20))
  ),

  tar_target(
    farm_temp_data_chunked, 
    command = {
      qs::qread(farm_ts_data_file) %>% 
        filter(farm_ID %in% farm_IDs_chunked[[1]]) %>%
        merge(farm_coords, by = "farm_ID") %>%
        mutate(keep = case_when(day >= t_start & day <= t_end ~ T, T ~ F)) %>% 
        dplyr::filter(keep) %>% 
        dplyr::select(farm_ID, day, temp_c)
    },
    pattern = farm_IDs_chunked
  ),

  tar_target(
    feed_names,
    command = {
      qs::qread(feed_params_file) %>% names()
    }
  ),
  tar_target(
    feed_params,
    qs::qread(feed_params_file)[[feed_names]],
    pattern = feed_names
  ),

  tar_target(
    reference_feed,
    qs::qread(feed_params_file)[["plant_dominant"]]
  ),

  tar_target(
    test_reference_feed,
    command = {
      fi <- sample(farm_IDs, 1)
      fg <- fish_growth(
        pop_params = all_params,
        species_params = all_params,
        water_temp = farm_temp_data_chunked$temp_c[farm_temp_data_chunked$farm_ID == fi],
        feed_params = reference_feed,
        times = c(
          t_start = min(farm_temp_data_chunked$day[farm_temp_data_chunked$farm_ID == fi]), 
          t_end = max(farm_temp_data_chunked$day[farm_temp_data_chunked$farm_ID == fi]), 
          dt = 1
        ),
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

  tar_target(all_params, c(qs::qread(species_params_file), qs::qread(pop_params_file))),

  # Calculate harvest size for each farm
  tar_target(
    harvest_size_chunked,
    command = {
      temp_data <- split(farm_temp_data_chunked, farm_temp_data_chunked$farm_ID)
      df <- purrr::map(temp_data, function(temp) {
        fg <- fish_growth(
          pop_params = all_params,
          species_params = all_params,
          water_temp = temp$temp_c,
          feed_params = reference_feed,
          times = c(t_start = min(temp$day), t_end = max(temp$day), dt = 1),
          init_weight = all_params["meanW"],
          ingmax = all_params["meanImax"]
        )
        fg %>%
          as.data.frame() %>%
          filter(!is.na(weight)) %>% 
          slice_tail(n = 1) %>%
          select(weight) %>% 
          mutate(farm_ID = unique(temp$farm_ID))
      })      
    },
    pattern = farm_temp_data_chunked
  ),

# Prepare parameters ----------------------------------------------------------------------------------------------
  tar_target(stat_names, levels(farm_run_chunked$measure)),

# Run growth ------------------------------------------------------------------------------------------------------
  tar_target(
    farm_static_data_chunked,
    command = {
      ts <- qs::qread(farm_ts_data_file) %>% 
        filter(farm_ID %in% farm_IDs_chunked[[1]]) %>% 
        merge(farm_coords, by = "farm_ID") %>%
        mutate(keep = case_when(day >= t_start & day <= t_end ~ T, T ~ F)) %>% 
        dplyr::filter(keep) %>% 
        dplyr::select(farm_ID, day, temp_c) %>% 
        group_by(farm_ID) %>% 
        reframe(
          t_start = min(day),
          t_end = max(day)
        )
      farm_production_file %>% 
        qs::qread() %>% 
        dplyr::select(farm_id, tonnes_per_farm) %>% 
        filter(farm_id %in% farm_IDs_chunked[[1]]) %>% 
        rename(farm_ID = farm_id) %>%
        distinct(farm_ID, tonnes_per_farm) %>% 
        mutate(prod = tonnes_per_farm %>% units::set_units("t") %>% units::set_units("g") %>% units::drop_units()) %>% 
        merge(ts, by = "farm_ID")
    },
    pattern = farm_IDs_chunked
  ),

  tar_target(
    farm_ts_data_chunked,
    command = {
      farm_temp_data <- split(farm_temp_data_chunked, farm_temp_data_chunked$farm_ID)
      farm_static_data <- split(farm_static_data_chunked, farm_static_data_chunked$farm_ID) %>% 
        purrr::map2(harvest_size_chunked, cbind)

      purrr::map2(farm_static_data, farm_temp_data, function(static, ts) {
        data.frame(
          Npop = generate_pop(
            harvest_n = static$prod/static$weight,
            mort = all_params['mortmyt'],
            times = c(t_start = static$t_start, t_end = static$t_end, dt = 1)
          ),
          temp_c = ts$temp_c
        )
      })
    },
    pattern = map(farm_static_data_chunked, harvest_size_chunked, farm_temp_data_chunked)
  ),


  tar_target(
    farm_run_chunked,
    command = {
      farm_static_data <- split(farm_static_data_chunked, farm_static_data_chunked$farm_ID)
      farm_ts_data <- farm_ts_data_chunked
      purrr::map2_dfr(farm_ts_data, farm_static_data, function(ts, static) {
        fg <- farm_growth(
          pop_params = all_params,
          species_params = all_params,
          water_temp = ts$temp_c,
          feed_params = feed_params,
          times = c(t_start = static$t_start, t_end = static$t_end, dt = 1),
          N_pop = ts$Npop,
          nruns = 5000
        )
        fg %>% 
          purrr::map_dfr(function(mat) {
            mat %>% 
              as.data.frame() %>% 
              rename(t = V1, mean = V2, sd = V3) %>% 
              mutate(
                farm_ID = as.integer(static$farm_ID),
                t = as.integer(t),
                feed = as.factor(feed_names)
              )
          }, .id = "measure") %>% 
          mutate(measure = as.factor(measure))
      })
    },
    pattern = cross(map(farm_ts_data_chunked,farm_static_data_chunked), map(feed_params, feed_names))
  ),

# Process growth ------------------------------------------------------------------------------------------------------
  tar_target(
    farm_results_chunked,
    command = {
      farm_run_chunked %>% 
        filter(measure == stat_names) %>% 
        mutate(prod_t = t-min(t)+1) %>% 
        filter(t != max(t))
    },
    pattern = cross(farm_run_chunked, stat_names)
  ),

  tar_target(
    biomass_produced_chunked,
    command = {
      farm_run_chunked %>% 
        filter(measure == "biomass_stat") %>% 
        mutate(prod_t = t-min(t)+1) %>% 
        filter(t == max(t))
    },
    pattern = farm_run_chunked
  ),

  tar_target(
    cohort_results_chunked,
    command = {
      farm_run <- split(farm_run_chunked, farm_run_chunked$farm_ID)
      purrr::map(farm_run, function(run) {
        cohort_1 <- run %>% filter(measure == stat_names & t != max(t)) # Fish are given no food on harvest day
        cohort_2 <- farm_to_cohort(cohort_1, time_offset = 365)
        cohort_3 <- farm_to_cohort(cohort_1, time_offset = 730)
        cohort_1 <- farm_to_cohort(cohort_1, time_offset = 0)
        lims <- min(cohort_2$t):(min(cohort_2$t)+365*2) # Two years after start of cohort 2

        tmp <- rbind(cohort_1, cohort_2, cohort_3) %>% 
          filter(t %in% lims) %>% 
          group_by(farm_ID, feed, t) %>%
          reframe(
            mean = sumna(mean),
            sd = sqrt(sumna(sd^2)),
            measure = as.factor(stat_names)
          ) %>% 
          mutate(prod_t = t-min(lims)+1)
      })
    },
    pattern = cross(farm_run_chunked, stat_names)
  )
)
