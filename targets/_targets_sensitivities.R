library(targets)
library(crew)
library(parallelly)

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "tibble"), 
  format = "qs", 
  controller = crew_controller_local(workers = 5, seconds_timeout = 90),
  workspace_on_error = TRUE
)

tar_source(c(
  "00_dirs.R",
  "00_model_functions.R"
))

list(
  # Load previously saved data --------------------------------------------------------------------------------
  tar_target(farm_coords_file, file.path(output_farm_data_path, "farm_coords.qs"), format = "file"),
  tar_target(farm_ts_data_file, file.path(output_farm_data_path, "farm_ts_data.qs"), format = "file"),
  tar_target(species_params_file, file.path(output_species_data_path, "species_params.qs"), format = "file"),
  tar_target(pop_params_file, file.path(output_species_data_path, "pop_params.qs"), format = "file"),
  tar_target(feed_params_file, file.path(output_species_data_path, "feed_params.qs"), format = "file"),

  tar_target(farm_coords, qs::qread(farm_coords_file)),
  tar_target(
    farm_IDs,
    qs::qread(farm_ts_data_file) %>% 
      distinct(farm_ID) %>% 
      pull(farm_ID)
  ),
  tar_target(
    farm_IDs_sample,
    farm_IDs %>% 
      sample(as.integer(length(farm_IDs)*0.1)) #%>% sample(3)
  ),
  tar_target(
    farm_ts_data, 
    qs::qread(farm_ts_data_file) %>% 
      merge(farm_coords, by = "farm_ID") %>% 
      dplyr::filter(farm_ID %in% farm_IDs) %>%
      mutate(keep = case_when(day >= t_start & day <= t_end ~ T, T ~ F)) %>% 
      dplyr::filter(keep) %>% 
      dplyr::select(farm_ID, day, temp_c)
  ),
  tar_target(
    reference_feed,
    qs::qread(feed_params_file)[["plant_dominant"]]
  ),

  # Get names of params to be adjusted
  tar_target(all_params, c(qs::qread(species_params_file), qs::qread(pop_params_file))),
  tar_target(param_names_pop, c("meanW", "deltaW", "meanImax", "deltaImax", "overFmean", "overFdelta")),
  tar_target(param_names_spec, names(all_params)[!names(all_params) %in% param_names_pop]),

  tar_target(factors, c(0.9, 1, 1.1)),
  tar_target(
    sens_run_spec, 
    command = {
      # Pre-split farm temperature data for speed
      farm_ts_list <- split(farm_ts_data, farm_ts_data$farm_ID)

      sens_farm <- purrr::map(factors, function(fact) {
        
        # Adjust params by factor
        adj_params <- all_params
        adj_params[param_names_spec] <- adj_params[param_names_spec] * fact

        purrr::map2(farm_IDs, farm_ts_list, function(farm, ts) {
          sens <- uni_farm_growth(
            pop_params = adj_params,
            species_params = adj_params,
            water_temp = ts$temp_c,
            feed_params = reference_feed,
            times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
            N_pop = rep(1, nrow(ts))
          )
          
          data.frame(
            weight = last(sens[['weight_stat']][, 2][!is.na(sens[['weight_stat']][, 2])]),
            dw = mean(sens[['dw_stat']][, 2], na.rm = T),
            total_excr = sum(sens[['total_excr_stat']][, 2], na.rm = T),
            P_excr = sum(sens[['P_excr_stat']][, 2], na.rm = T),
            L_excr = sum(sens[['L_excr_stat']][, 2], na.rm = T),
            C_excr = sum(sens[['C_excr_stat']][, 2], na.rm = T),
            total_uneat = sum(sens[['total_uneat_stat']][, 2], na.rm = T),
            P_uneat = sum(sens[['P_uneat_stat']][, 2], na.rm = T),
            L_uneat = sum(sens[['L_uneat_stat']][, 2], na.rm = T),
            C_uneat = sum(sens[['C_uneat_stat']][, 2], na.rm = T),
            food_prov = sum(sens[['food_prov_stat']][, 2], na.rm = T),
            rel_feeding = mean(sens[['rel_feeding_stat']][, 2], na.rm = T),
            ing_act = sum(sens[['ing_act_stat']][, 2], na.rm = T),
            anab = mean(sens[['anab_stat']][, 2], na.rm = T),
            catab = mean(sens[['catab_stat']][, 2], na.rm = T),
            O2 = sum(sens[['O2_stat']][, 2], na.rm = T),
            NH4 = sum(sens[['NH4_stat']][, 2], na.rm = T),
            adj_param = as.factor(param_names_spec),
            farm_ID = farm
          )
        }) %>% bind_rows() %>% 
          mutate(factor = fact)
      }) %>% bind_rows()
    },
    pattern = param_names_spec
  ),
  
  tar_target(
    sens_results_spec,
    command = {
      cols <- colnames(sens_run_spec)[!colnames(sens_run_spec) %in% c("adj_param", "farm_ID", "factor")]
      sens <- sens_run_spec %>%
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
    pattern = sens_run_spec
  ),
  
  tar_target(
    sens_run_pop, 
    command = {
      # Pre-split farm temperature data for speed
      farm_ts_list <- filter(farm_ts_data, farm_ID %in% farm_IDs_sample)
      farm_ts_list <- split(farm_ts_list, farm_ts_list$farm_ID)

      purrr::map(factors, function(fact) {
        # Adjust params by factor
        adj_params <- all_params
        adj_params[param_names_pop] <- adj_params[param_names_pop] * fact

        purrr::map2(farm_IDs_sample, farm_ts_list, function(farm, ts) {
          sens <- farm_growth(
            pop_params = adj_params,
            species_params = adj_params,
            water_temp = ts$temp_c,
            feed_params = reference_feed,
            times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
            N_pop = rep(1, nrow(ts)),
            nruns = 500
          )
          
          data.frame(
            weight = last(sens[['weight_stat']][, 2][!is.na(sens[['weight_stat']][, 2])]),
            dw = mean(sens[['dw_stat']][, 2], na.rm = T),
            total_excr = sum(sens[['total_excr_stat']][, 2], na.rm = T),
            P_excr = sum(sens[['P_excr_stat']][, 2], na.rm = T),
            L_excr = sum(sens[['L_excr_stat']][, 2], na.rm = T),
            C_excr = sum(sens[['C_excr_stat']][, 2], na.rm = T),
            total_uneat = sum(sens[['total_uneat_stat']][, 2], na.rm = T),
            P_uneat = sum(sens[['P_uneat_stat']][, 2], na.rm = T),
            L_uneat = sum(sens[['L_uneat_stat']][, 2], na.rm = T),
            C_uneat = sum(sens[['C_uneat_stat']][, 2], na.rm = T),
            food_prov = sum(sens[['food_prov_stat']][, 2], na.rm = T),
            rel_feeding = mean(sens[['rel_feeding_stat']][, 2], na.rm = T),
            ing_act = sum(sens[['ing_act_stat']][, 2], na.rm = T),
            anab = mean(sens[['anab_stat']][, 2], na.rm = T),
            catab = mean(sens[['catab_stat']][, 2], na.rm = T),
            O2 = sum(sens[['O2_stat']][, 2], na.rm = T),
            NH4 = sum(sens[['NH4_stat']][, 2], na.rm = T),
            adj_param = as.factor(param_names_pop),
            farm_ID = as.integer(farm)
          )
        }) %>% 
          bind_rows() %>% 
          mutate(factor = fact)
      }) %>% 
        bind_rows()
    },
    pattern = param_names_pop
  ),

  tar_target(
    sens_results_pop,
    command = {
      cols <- colnames(sens_run_pop)[!colnames(sens_run_pop) %in% c("adj_param", "farm_ID", "factor")]
      sens <- sens_run_pop %>%
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
    pattern = sens_run_pop
  )
)
