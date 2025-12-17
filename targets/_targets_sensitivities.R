suppressPackageStartupMessages(suppressWarnings({
  library(targets)
  library(crew)
  library(magrittr)
  library(qs)
  library(tidyverse)
  # library(arrow)
  # library(tibble)
}))

tar_option_set(
  format = "qs", 
  controller = crew::crew_controller_local(
    workers = parallelly::availableCores()-1, 
    seconds_idle = 120
  ),
  workspace_on_error = TRUE,
  garbage_collection = 10
)

tar_source(
  files = c(
    "/home/treimer/local_to_global_mariculture_modelling/src/dirs.R",
    "/home/treimer/local_to_global_mariculture_modelling/src/functions.R",
    "/home/treimer/local_to_global_mariculture_modelling/src/model_functions.R"
  )
)

# Dirs
output_path <- here() %>% file.path("outputs")
output_species_data_path <- file.path(output_path, "species_data")

# Globals
inds_per_farm <- 250
farm_sample <- 272
farm_chunk_size <- 25
reference_feed_name <- "marine_dominant_biomar"

list(
  # Load previously saved data --------------------------------------------------------------------------------
  tar_target(farm_coords_file, file.path(output_farm_data_path, "farm_coords.qs"), format = "file"),
  tar_target(farm_ts_data_file, file.path(output_farm_data_path, "farm_ts_data.qs"), format = "file"),
  tar_target(species_params_file, file.path(output_species_data_path, "species_params.qs"), format = "file"),
  tar_target(pop_params_file, file.path(output_species_data_path, "pop_params.qs"), format = "file"),
  tar_target(feed_params_file, file.path(output_species_data_path, "feed_params.qs"), format = "file"),

  tar_target(farm_coords, qread(farm_coords_file)),
  tar_target(
    farm_IDs,
    command = {
      qread(farm_ts_data_file) %>% 
        filter(farm_ID != 2703) %>% # The only Russian farm?
        distinct(farm_ID) %>% 
        pull(farm_ID) %>% 
        sample(farm_sample)
    }
  ),

  tar_target(
    farm_IDs_chunked,
    split(farm_IDs, ceiling(seq_along(farm_IDs)/farm_chunk_size))
  ),

  tar_target(
    farm_ts_data, 
    command = {
      qread(farm_ts_data_file) %>% 
        merge(farm_coords, by = "farm_ID") %>% 
        dplyr::filter(farm_ID %in% farm_IDs) %>%
        mutate(keep = case_when(day >= t_start & day <= t_end ~ T, T ~ F)) %>% 
        dplyr::filter(keep) %>% 
        dplyr::select(farm_ID, day, temp_c)
    }
  ),

  tar_target(
    farm_temp_data_chunked, 
    farm_ts_data %>% filter(farm_ID %in% farm_IDs_chunked[[1]]),
    pattern = farm_IDs_chunked
  ),

    tar_target(
    test_reference_feed,
    command = {
      fi <- sample(farm_IDs, 1)
      fg <- fish_growth(
        pop_params = all_params,
        species_params = all_params,
        water_temp = farm_ts_data$temp_c[farm_ts_data$farm_ID == fi],
        feed_params = qread(feed_params_file)[[reference_feed_name]],
        times = c(
          t_start = min(farm_ts_data$day[farm_ts_data$farm_ID == fi]), 
          t_end = max(farm_ts_data$day[farm_ts_data$farm_ID == fi]), 
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

  # Get names of params to be adjusted
  tar_target(all_params, c(qread(species_params_file), qread(pop_params_file))),
  tar_target(param_names_pop, c("meanW", "deltaW", "meanImax", "deltaImax", "overFmean", "overFdelta")),
  tar_target(param_names_spec, names(all_params)[!names(all_params) %in% param_names_pop]),

  tar_target(factors, c(0.9, 1, 1.1)),
  tar_target(
    sens_run_spec_lo, 
    command = {
      # Pre-split farm temperature data for speed
      farm_ts_list <- split(farm_ts_data, farm_ts_data$farm_ID)
        
      # Adjust params by factor
      adj_params <- all_params
      adj_params[param_names_spec] <- adj_params[param_names_spec] * factors[1]

      map2_dfr(farm_IDs, farm_ts_list, function(farm, ts) {
        fg <- uni_farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = qread(feed_params_file)[[reference_feed_name]],
          times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(ts))
        )
        
        data.frame(
          weight = last(fg[['weight_stat']][, 2][!is.na(fg[['weight_stat']][, 2])]),
          dw = mean(fg[['dw_stat']][, 2], na.rm = T),
          total_excr = sum(fg[['total_excr_stat']][, 2], na.rm = T),
          P_excr = sum(fg[['P_excr_stat']][, 2], na.rm = T),
          L_excr = sum(fg[['L_excr_stat']][, 2], na.rm = T),
          C_excr = sum(fg[['C_excr_stat']][, 2], na.rm = T),
          total_uneat = sum(fg[['total_uneat_stat']][, 2], na.rm = T),
          P_uneat = sum(fg[['P_uneat_stat']][, 2], na.rm = T),
          L_uneat = sum(fg[['L_uneat_stat']][, 2], na.rm = T),
          C_uneat = sum(fg[['C_uneat_stat']][, 2], na.rm = T),
          food_prov = sum(fg[['food_prov_stat']][, 2], na.rm = T),
          rel_feeding = mean(fg[['rel_feeding_stat']][, 2], na.rm = T),
          ing_act = sum(fg[['ing_act_stat']][, 2], na.rm = T),
          anab = mean(fg[['anab_stat']][, 2], na.rm = T),
          catab = mean(fg[['catab_stat']][, 2], na.rm = T),
          O2 = sum(fg[['O2_stat']][, 2], na.rm = T),
          NH4 = sum(fg[['NH4_stat']][, 2], na.rm = T),
          adj_param = as.factor(param_names_spec),
          farm_ID = farm,
          factor = factors[1]
        )
      })
    },
    pattern = param_names_spec
  ),

  tar_target(
    sens_run_spec_mi, 
    command = {
      # Pre-split farm temperature data for speed
      farm_ts_list <- split(farm_ts_data, farm_ts_data$farm_ID)
        
      # Adjust params by factor
      adj_params <- all_params
      adj_params[param_names_spec] <- adj_params[param_names_spec] * factors[2]

      map2_dfr(farm_IDs, farm_ts_list, function(farm, ts) {
        fg <- uni_farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = qread(feed_params_file)[[reference_feed_name]],
          times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(ts))
        )
        
        data.frame(
          weight = last(fg[['weight_stat']][, 2][!is.na(fg[['weight_stat']][, 2])]),
          dw = mean(fg[['dw_stat']][, 2], na.rm = T),
          total_excr = sum(fg[['total_excr_stat']][, 2], na.rm = T),
          P_excr = sum(fg[['P_excr_stat']][, 2], na.rm = T),
          L_excr = sum(fg[['L_excr_stat']][, 2], na.rm = T),
          C_excr = sum(fg[['C_excr_stat']][, 2], na.rm = T),
          total_uneat = sum(fg[['total_uneat_stat']][, 2], na.rm = T),
          P_uneat = sum(fg[['P_uneat_stat']][, 2], na.rm = T),
          L_uneat = sum(fg[['L_uneat_stat']][, 2], na.rm = T),
          C_uneat = sum(fg[['C_uneat_stat']][, 2], na.rm = T),
          food_prov = sum(fg[['food_prov_stat']][, 2], na.rm = T),
          rel_feeding = mean(fg[['rel_feeding_stat']][, 2], na.rm = T),
          ing_act = sum(fg[['ing_act_stat']][, 2], na.rm = T),
          anab = mean(fg[['anab_stat']][, 2], na.rm = T),
          catab = mean(fg[['catab_stat']][, 2], na.rm = T),
          O2 = sum(fg[['O2_stat']][, 2], na.rm = T),
          NH4 = sum(fg[['NH4_stat']][, 2], na.rm = T),
          adj_param = as.factor(param_names_spec),
          farm_ID = farm,
          factor = factors[2]
        )
      })
    },
    pattern = param_names_spec
  ),

    tar_target(
    sens_run_spec_hi, 
    command = {
      # Pre-split farm temperature data for speed
      farm_ts_list <- split(farm_ts_data, farm_ts_data$farm_ID)
        
      # Adjust params by factor
      adj_params <- all_params
      adj_params[param_names_spec] <- adj_params[param_names_spec] * factors[3]

      map2_dfr(farm_IDs, farm_ts_list, function(farm, ts) {
        fg <- uni_farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = qread(feed_params_file)[[reference_feed_name]],
          times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(ts))
        )
        
        data.frame(
          weight = last(fg[['weight_stat']][, 2][!is.na(fg[['weight_stat']][, 2])]),
          dw = mean(fg[['dw_stat']][, 2], na.rm = T),
          total_excr = sum(fg[['total_excr_stat']][, 2], na.rm = T),
          P_excr = sum(fg[['P_excr_stat']][, 2], na.rm = T),
          L_excr = sum(fg[['L_excr_stat']][, 2], na.rm = T),
          C_excr = sum(fg[['C_excr_stat']][, 2], na.rm = T),
          total_uneat = sum(fg[['total_uneat_stat']][, 2], na.rm = T),
          P_uneat = sum(fg[['P_uneat_stat']][, 2], na.rm = T),
          L_uneat = sum(fg[['L_uneat_stat']][, 2], na.rm = T),
          C_uneat = sum(fg[['C_uneat_stat']][, 2], na.rm = T),
          food_prov = sum(fg[['food_prov_stat']][, 2], na.rm = T),
          rel_feeding = mean(fg[['rel_feeding_stat']][, 2], na.rm = T),
          ing_act = sum(fg[['ing_act_stat']][, 2], na.rm = T),
          anab = mean(fg[['anab_stat']][, 2], na.rm = T),
          catab = mean(fg[['catab_stat']][, 2], na.rm = T),
          O2 = sum(fg[['O2_stat']][, 2], na.rm = T),
          NH4 = sum(fg[['NH4_stat']][, 2], na.rm = T),
          adj_param = as.factor(param_names_spec),
          farm_ID = farm,
          factor = factors[3]
        )
      })
    },
    pattern = param_names_spec
  ),

  tar_target(
    sens_results_spec,
    command = {
      sens <- rbind(
        sens_run_spec_lo,
        sens_run_spec_mi,
        sens_run_spec_hi
      ) 
      cols <- colnames(sens)[!colnames(sens) %in% c("adj_param", "farm_ID", "factor")]
      
      sens <- sens %>%
        pivot_longer(
          cols = all_of(cols), 
          names_to = "measure", 
          values_to = "value", 
          names_transform = list(measure = as.factor)
        ) %>% 
        pivot_wider(
          names_from = "factor", 
          values_from = "value", 
          names_prefix = "fact_"
        ) %>% 
        mutate(sensitivity = (fact_1.1-fact_0.9)/(0.2*fact_1)) 
      
      sens %>% 
        group_by(adj_param, measure) %>%
        reframe(
          mean_sens = meanna(sensitivity),
          sd_sens = sdna(sensitivity)
        )
    },
    pattern = map(sens_run_spec_lo, sens_run_spec_hi, sens_run_spec_mi)
  ),
  
  tar_target(
    sens_run_pop_lo, 
    command = {
      # Adjust params by factor
      adj_params <- all_params
      adj_params[param_names_pop] <- adj_params[param_names_pop] * factors[1]

      farm_ts <- split(farm_temp_data_chunked, farm_temp_data_chunked$farm_ID)
      # tic()
      map2_dfr(names(farm_ts), farm_ts, function(farm_ID, ts) {
        fg <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = qread(feed_params_file)[[reference_feed_name]],
          times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(ts)),
          nruns = inds_per_farm
        )
      
        data.frame(
          weight = last(fg[['weight_stat']][, 2][!is.na(fg[['weight_stat']][, 2])]),
          dw = mean(fg[['dw_stat']][, 2], na.rm = T),
          total_excr = sum(fg[['total_excr_stat']][, 2], na.rm = T),
          P_excr = sum(fg[['P_excr_stat']][, 2], na.rm = T),
          L_excr = sum(fg[['L_excr_stat']][, 2], na.rm = T),
          C_excr = sum(fg[['C_excr_stat']][, 2], na.rm = T),
          total_uneat = sum(fg[['total_uneat_stat']][, 2], na.rm = T),
          P_uneat = sum(fg[['P_uneat_stat']][, 2], na.rm = T),
          L_uneat = sum(fg[['L_uneat_stat']][, 2], na.rm = T),
          C_uneat = sum(fg[['C_uneat_stat']][, 2], na.rm = T),
          food_prov = sum(fg[['food_prov_stat']][, 2], na.rm = T),
          rel_feeding = mean(fg[['rel_feeding_stat']][, 2], na.rm = T),
          ing_act = sum(fg[['ing_act_stat']][, 2], na.rm = T),
          anab = mean(fg[['anab_stat']][, 2], na.rm = T),
          catab = mean(fg[['catab_stat']][, 2], na.rm = T),
          O2 = sum(fg[['O2_stat']][, 2], na.rm = T),
          NH4 = sum(fg[['NH4_stat']][, 2], na.rm = T),
          adj_param = as.factor(param_names_pop),
          farm_ID = as.integer(farm_ID),
          factor = factors[1]
        )
      })
      # toc()
    },
    pattern = cross(farm_temp_data_chunked, param_names_pop)
  ),

  tar_target(
    sens_run_pop_mi, 
    command = {
      # Adjust params by factor
      adj_params <- all_params
      adj_params[param_names_pop] <- adj_params[param_names_pop] * factors[2]

      farm_ts <- split(farm_temp_data_chunked, farm_temp_data_chunked$farm_ID)
      map2_dfr(names(farm_ts), farm_ts, function(farm_ID, ts) {
        fg <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = qread(feed_params_file)[[reference_feed_name]],
          times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(ts)),
          nruns = inds_per_farm
        )
      
        data.frame(
          weight = last(fg[['weight_stat']][, 2][!is.na(fg[['weight_stat']][, 2])]),
          dw = mean(fg[['dw_stat']][, 2], na.rm = T),
          total_excr = sum(fg[['total_excr_stat']][, 2], na.rm = T),
          P_excr = sum(fg[['P_excr_stat']][, 2], na.rm = T),
          L_excr = sum(fg[['L_excr_stat']][, 2], na.rm = T),
          C_excr = sum(fg[['C_excr_stat']][, 2], na.rm = T),
          total_uneat = sum(fg[['total_uneat_stat']][, 2], na.rm = T),
          P_uneat = sum(fg[['P_uneat_stat']][, 2], na.rm = T),
          L_uneat = sum(fg[['L_uneat_stat']][, 2], na.rm = T),
          C_uneat = sum(fg[['C_uneat_stat']][, 2], na.rm = T),
          food_prov = sum(fg[['food_prov_stat']][, 2], na.rm = T),
          rel_feeding = mean(fg[['rel_feeding_stat']][, 2], na.rm = T),
          ing_act = sum(fg[['ing_act_stat']][, 2], na.rm = T),
          anab = mean(fg[['anab_stat']][, 2], na.rm = T),
          catab = mean(fg[['catab_stat']][, 2], na.rm = T),
          O2 = sum(fg[['O2_stat']][, 2], na.rm = T),
          NH4 = sum(fg[['NH4_stat']][, 2], na.rm = T),
          adj_param = as.factor(param_names_pop),
          farm_ID = as.integer(farm_ID),
          factor = factors[2]
        )
      })
    },
    pattern = cross(farm_temp_data_chunked, param_names_pop)
  ),

  tar_target(
    sens_run_pop_hi, 
    command = {
      # Adjust params by factor
      adj_params <- all_params
      adj_params[param_names_pop] <- adj_params[param_names_pop] * factors[3]

      farm_ts <- split(farm_temp_data_chunked, farm_temp_data_chunked$farm_ID)
      map2_dfr(names(farm_ts), farm_ts, function(farm_ID, ts) {
        fg <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = qread(feed_params_file)[[reference_feed_name]],
          times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(ts)),
          nruns = inds_per_farm
        )
      
        data.frame(
          weight = last(fg[['weight_stat']][, 2][!is.na(fg[['weight_stat']][, 2])]),
          dw = mean(fg[['dw_stat']][, 2], na.rm = T),
          total_excr = sum(fg[['total_excr_stat']][, 2], na.rm = T),
          P_excr = sum(fg[['P_excr_stat']][, 2], na.rm = T),
          L_excr = sum(fg[['L_excr_stat']][, 2], na.rm = T),
          C_excr = sum(fg[['C_excr_stat']][, 2], na.rm = T),
          total_uneat = sum(fg[['total_uneat_stat']][, 2], na.rm = T),
          P_uneat = sum(fg[['P_uneat_stat']][, 2], na.rm = T),
          L_uneat = sum(fg[['L_uneat_stat']][, 2], na.rm = T),
          C_uneat = sum(fg[['C_uneat_stat']][, 2], na.rm = T),
          food_prov = sum(fg[['food_prov_stat']][, 2], na.rm = T),
          rel_feeding = mean(fg[['rel_feeding_stat']][, 2], na.rm = T),
          ing_act = sum(fg[['ing_act_stat']][, 2], na.rm = T),
          anab = mean(fg[['anab_stat']][, 2], na.rm = T),
          catab = mean(fg[['catab_stat']][, 2], na.rm = T),
          O2 = sum(fg[['O2_stat']][, 2], na.rm = T),
          NH4 = sum(fg[['NH4_stat']][, 2], na.rm = T),
          adj_param = as.factor(param_names_pop),
          farm_ID = as.integer(farm_ID),
          factor = factors[3]
        )
      })
    },
    pattern = cross(farm_temp_data_chunked, param_names_pop)
  ),

  tar_target(
    sens_results_pop_0,
    command = {
      sens <- rbind(sens_run_pop_lo, sens_run_pop_mi, sens_run_pop_hi) 
      cols <- colnames(sens)[!colnames(sens) %in% c("adj_param", "farm_ID", "factor")]
      
      sens <- sens %>%
        pivot_longer(
          cols = all_of(cols), 
          names_to = "measure", 
          values_to = "value", 
          names_transform = list(measure = as.factor)
        ) %>% 
        pivot_wider(
          names_from = "factor", 
          values_from = "value", 
          names_prefix = "fact_"
        ) %>% 
        mutate(sensitivity = (fact_1.1-fact_0.9)/(0.2*fact_1))
    },
    pattern = map(sens_run_pop_lo, sens_run_pop_mi, sens_run_pop_hi)
  ),
  
  tar_target(
    sens_results_pop,
    command = {
      sens_results_pop_0 %>% 
        group_by(adj_param, measure) %>%
        reframe(
          mean_sens = meanna(sensitivity),
          sd_sens = sdna(sensitivity)
        )
    }
  )
)
