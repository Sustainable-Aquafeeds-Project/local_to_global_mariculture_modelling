library(targets)
library(crew)
library(parallelly)

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "tibble"), 
  format = "qs", 
  controller = crew_controller_local(workers = 18, seconds_timeout = 90),
  workspace_on_error = TRUE
)

list(
  tar_target(
    tar_farm_IDs,
    farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID)
  ),
  tar_target(
    tar_farm_ID_sample,
    farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID) %>% sample(25) # %>% sample(271)
  ),
  
  tar_target(tar_param_names_pop, c("meanW", "deltaW", "meanImax", "deltaImax", "overFmean", "overFdelta")),
  tar_target(tar_param_names_spec, sens_params_names[!sens_params_names %in% tar_param_names_pop]),
  
  tar_target(
    tar_farm_ts,
    command = {
      t_start <- farm_coords$t_start[farm_coords$farm_ID == tar_farm_IDs]
      t_end <- farm_coords$t_end[farm_coords$farm_ID == tar_farm_IDs]
      farm_ts_data %>%
        dplyr::filter(farm_ID == tar_farm_IDs) %>%
        dplyr::filter(between(day, t_start, t_end))
    },
    pattern = tar_farm_IDs
  ),
  
  tar_target(
    tar_sens_run_spec_lo, 
    command = {
      tar_farm_ts_list <- split(tar_farm_ts, tar_farm_ts$farm_ID)
      sens_farm <- purrr::map2_dfr(tar_farm_IDs, tar_farm_ts_list, function(farm, ts) {
        adj_params <- sens_all_params
        adj_params[tar_param_names_spec] <- adj_params[tar_param_names_spec] * factors[1]
        
        sens <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = reference_feed,
          times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(ts)),
          nruns = 10
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
          adj_param = as.factor(tar_param_names_spec)
        ) %>%
          mutate(farm_ID = farm)
      }) %>%
        mutate(factor = factors[1])
      
      # colnms <- colnames(sens_farm)[!colnames(sens_farm) %in% c("adj_param", "factor", "farm_ID")]
      # df <- sens_farm %>% 
      #   pivot_longer(cols = all_of(colnms), values_to = "value", 
      #                names_to = "measure", names_transform = list(measure = as.factor)) %>%
      #   group_by(adj_param, factor, farm_ID, measure) %>%
      #   pivot_wider(names_from = "factor", names_prefix = "factor_", values_from = "value") %>%
      #   mutate(sensitivity = (factor_1.1 - factor_0.9) / (0.2 * factor_1)) %>%
      #   group_by(adj_param, measure, farm_ID) %>%
      #   reframe(sensitivity = meanna(sensitivity))
    },
    pattern = tar_param_names_spec
  ),
  
  tar_target(
    tar_sens_run_spec_mi, 
    command = {
      tar_farm_ts_list <- split(tar_farm_ts, tar_farm_ts$farm_ID)
      sens_farm <- purrr::map2_dfr(tar_farm_IDs, tar_farm_ts_list, function(farm, ts) {
        adj_params <- sens_all_params
        adj_params[tar_param_names_spec] <- adj_params[tar_param_names_spec] * factors[2]
        
        sens <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = reference_feed,
          times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(ts)),
          nruns = 10
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
          adj_param = as.factor(tar_param_names_spec)
        ) %>%
          mutate(farm_ID = farm)
      }) %>%
        mutate(factor = factors[2])
    },
    pattern = tar_param_names_spec
  ),
  
  tar_target(
    tar_sens_run_spec_hi, 
    command = {
      tar_farm_ts_list <- split(tar_farm_ts, tar_farm_ts$farm_ID)
      sens_farm <- purrr::map2_dfr(tar_farm_IDs, tar_farm_ts_list, function(farm, ts) {
        adj_params <- sens_all_params
        adj_params[tar_param_names_spec] <- adj_params[tar_param_names_spec] * factors[3]
        
        sens <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = ts$temp_c,
          feed_params = reference_feed,
          times = c('t_start' = min(ts$day), 't_end' = max(ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(ts)),
          nruns = 10
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
          adj_param = as.factor(tar_param_names_spec)
        ) %>%
          mutate(farm_ID = farm)
      }) %>%
        mutate(factor = factors[3])
    },
    pattern = tar_param_names_spec
  ),
  
  tar_target(
    tar_sens_results_spec,
    command = {
      cols <- colnames(tar_sens_run_spec_lo)[
        !colnames(tar_sens_run_spec_lo) %in% c("adj_param", "farm_ID", "factor")
      ]
      rbind(
        tar_sens_run_spec_lo,
        tar_sens_run_spec_mi,
        tar_sens_run_spec_hi
      ) %>%
        pivot_longer(cols = all_of(cols), names_to = "measure", values_to = "value", 
                     names_transform = list(measure = as.factor)) %>% 
        pivot_wider(names_from = "factor", values_from = "value", names_prefix = "fact_") %>% 
        mutate(sensitivity = (fact_1.1-fact_0.9)/(0.2*fact_1)) %>% 
        group_by(adj_param, measure) %>%
        reframe(mean_sens = meanna(sensitivity),
                sd_sens = sdna(sensitivity))
    },
    pattern = map(tar_sens_run_spec_lo, tar_sens_run_spec_mi, tar_sens_run_spec_hi)
  ),
  
  tar_target(
    tar_sens_run_pop_lo, 
    command = {
      adj_params <- sens_all_params
      adj_params[tar_param_names_pop] <- adj_params[tar_param_names_pop] * factors[1]
      
      ts_temp <- tar_farm_ts[tar_farm_ts$farm_ID == tar_farm_ID_sample, ]
      
      sens <- farm_growth(
        pop_params = adj_params,
        species_params = adj_params,
        water_temp = ts_temp$temp_c,
        feed_params = reference_feed,
        times = c('t_start' = min(ts_temp$day), 't_end' = max(ts_temp$day), 'dt' = 1),
        N_pop = rep(1, nrow(ts_temp)),
        nruns = 5000
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
        adj_param = as.factor(tar_param_names_pop)
      ) %>%
        mutate(farm_ID = tar_farm_ID_sample,
               factor = factors[1])
    },
    pattern = cross(tar_param_names_pop, tar_farm_ID_sample)
  ),

  tar_target(
    tar_sens_run_pop_mi, 
    command = {
      adj_params <- sens_all_params
      adj_params[tar_param_names_pop] <- adj_params[tar_param_names_pop] * factors[2]
      
      ts_temp <- tar_farm_ts[tar_farm_ts$farm_ID == tar_farm_ID_sample, ]
      
      sens <- farm_growth(
        pop_params = adj_params,
        species_params = adj_params,
        water_temp = ts_temp$temp_c,
        feed_params = reference_feed,
        times = c('t_start' = min(ts_temp$day), 't_end' = max(ts_temp$day), 'dt' = 1),
        N_pop = rep(1, nrow(ts_temp)),
        nruns = 5000
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
        adj_param = as.factor(tar_param_names_pop)
      ) %>%
        mutate(farm_ID = tar_farm_ID_sample,
               factor = factors[2])
    },
    pattern = cross(tar_param_names_pop, tar_farm_ID_sample)
  ),

  tar_target(
    tar_sens_run_pop_hi, 
    command = {
      adj_params <- sens_all_params
      adj_params[tar_param_names_pop] <- adj_params[tar_param_names_pop] * factors[3]
      
      ts_temp <- tar_farm_ts[tar_farm_ts$farm_ID == tar_farm_ID_sample, ]
      
      sens <- farm_growth(
        pop_params = adj_params,
        species_params = adj_params,
        water_temp = ts_temp$temp_c,
        feed_params = reference_feed,
        times = c('t_start' = min(ts_temp$day), 't_end' = max(ts_temp$day), 'dt' = 1),
        N_pop = rep(1, nrow(ts_temp)),
        nruns = 5000
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
        adj_param = as.factor(tar_param_names_pop)
      ) %>%
        mutate(farm_ID = tar_farm_ID_sample,
               factor = factors[3])
    },
    pattern = cross(tar_param_names_pop, tar_farm_ID_sample)
  ),
  
  tar_target(
    tar_sens_results_pop,
    command = {
      cols <- colnames(tar_sens_run_pop_lo)[
        !colnames(tar_sens_run_pop_lo) %in% c("adj_param", "farm_ID", "factor")
      ]
      rbind(
        tar_sens_run_pop_lo,
        tar_sens_run_pop_mi,
        tar_sens_run_pop_hi
      ) %>%
        pivot_longer(cols = all_of(cols), names_to = "measure", values_to = "value", 
                     names_transform = list(measure = as.factor)) %>% 
        pivot_wider(names_from = "factor", values_from = "value", names_prefix = "fact_") %>% 
        mutate(sensitivity = (fact_1.1-fact_0.9)/(0.2*fact_1)) %>% 
        group_by(adj_param, measure) %>%
        reframe(mean_sens = meanna(sensitivity),
                sd_sens = sdna(sensitivity))
    },
    pattern = map(tar_sens_run_pop_lo, tar_sens_run_pop_mi, tar_sens_run_pop_hi)
  )
)
