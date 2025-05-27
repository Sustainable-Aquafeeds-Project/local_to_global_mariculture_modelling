library(targets)
library(crew)
library(parallelly)

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "tibble"), 
  format = "qs", 
  controller = crew_controller_local(workers = 10),
  workspace_on_error = TRUE
)

list(
  tar_target(
    tar_farm_IDs,
    farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID)
  ),
  
  tar_target(tar_factors, factors),
  tar_target(tar_param_names_pop, c("meanW", "deltaW", "meanImax", "deltaImax", "overFmean", "overFdelta")),
  tar_target(
    tar_param_names_spec,
    sens_params_names[!sens_params_names %in% tar_param_names_pop] %>% 
    sample(4)
  ),
  
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
    tar_sens_run_spec, 
    command = {
      purrr::map_dfr(factors, function(fact) {
        adj_params <- sens_all_params
        adj_params[tar_param_names_spec] <- adj_params[tar_param_names_spec] * fact
        
        sens <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = tar_farm_ts$temp_c,
          feed_params = reference_feed,
          times = c('t_start' = min(tar_farm_ts$day), 't_end' = max(tar_farm_ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(tar_farm_ts)),
          nruns = 100
        )
      
        data.frame(
          weight = last(sens[['weight_stat']][,2][!is.na(sens[['weight_stat']][,2])]),
          dw = mean(sens[['dw_stat']][,2], na.rm = T),
          total_excr = sum(sens[['total_excr_stat']][,2], na.rm = T),
          P_excr = sum(sens[['P_excr_stat']][,2], na.rm = T),
          L_excr = sum(sens[['L_excr_stat']][,2], na.rm = T),
          C_excr = sum(sens[['C_excr_stat']][,2], na.rm = T),
          total_uneat = sum(sens[['total_uneat_stat']][,2], na.rm = T),
          P_uneat = sum(sens[['P_uneat_stat']][,2], na.rm = T),
          L_uneat = sum(sens[['L_uneat_stat']][,2], na.rm = T),
          C_uneat = sum(sens[['C_uneat_stat']][,2], na.rm = T),
          food_prov = sum(sens[['food_prov_stat']][,2], na.rm = T),
          rel_feeding = mean(sens[['rel_feeding_stat']][,2], na.rm = T),
          ing_act = sum(sens[['ing_act_stat']][,2], na.rm = T),
          anab = mean(sens[['anab_stat']][,2], na.rm = T),
          catab = mean(sens[['catab_stat']][,2], na.rm = T),
          O2 = sum(sens[['O2_stat']][,2], na.rm = T),
          NH4 = sum(sens[['NH4_stat']][,2], na.rm = T),
          adj_param = as.factor(tar_param_names_spec),
          factor = fact,
          farm_ID = tar_farm_IDs
        )
      }) %>%
        pivot_longer(
          cols = c(weight, dw, total_excr, P_excr, L_excr, C_excr, total_uneat, P_uneat, L_uneat, C_uneat,
                   food_prov, rel_feeding, ing_act, anab, catab, O2, NH4),
          names_to = "measure", names_transform = list(measure = as.factor), values_to = "value"
        ) %>%
        group_by(adj_param, factor, farm_ID, measure) %>%
        pivot_wider(names_from = "factor", names_prefix = "fact_", values_from = "value") %>%
        mutate(sensitivity = (fact_1.1 - fact_0.9) / (0.2 * fact_1)) %>%
        group_by(adj_param, measure, farm_ID) %>%
        reframe(sensitivity = meanna(sensitivity))
    },
    pattern = cross(tar_param_names_spec, map(tar_farm_ts, tar_farm_IDs))
  ),
  
  tar_target(
    tar_sens_results_spec,
    command = {
      tar_sens_run_spec %>%
        group_by(adj_param, measure) %>%
        reframe(mean_sens = meanna(sensitivity),
                sd_sens = sdna(sensitivity))
    },
    pattern = tar_param_names_spec
  ),
  
  tar_target(
    tar_sens_run_pop, 
    command = {
      purrr::map_dfr(factors, function(fact) {
        adj_params <- sens_all_params
        adj_params[tar_param_names_pop] <- adj_params[tar_param_names_pop] * factors
        
        sens <- farm_growth(
          pop_params = adj_params,
          species_params = adj_params,
          water_temp = tar_farm_ts$temp_c,
          feed_params = reference_feed,
          times = c('t_start' = min(tar_farm_ts$day), 't_end' = max(tar_farm_ts$day), 'dt' = 1),
          N_pop = rep(1, nrow(tar_farm_ts)),
          nruns = 5000
        )
        
        data.frame(
          weight = last(sens[['weight_stat']][,2][!is.na(sens[['weight_stat']][,2])]),
          dw = mean(sens[['dw_stat']][,2], na.rm = T),
          total_excr = sum(sens[['total_excr_stat']][,2], na.rm = T),
          P_excr = sum(sens[['P_excr_stat']][,2], na.rm = T),
          L_excr = sum(sens[['L_excr_stat']][,2], na.rm = T),
          C_excr = sum(sens[['C_excr_stat']][,2], na.rm = T),
          total_uneat = sum(sens[['total_uneat_stat']][,2], na.rm = T),
          P_uneat = sum(sens[['P_uneat_stat']][,2], na.rm = T),
          L_uneat = sum(sens[['L_uneat_stat']][,2], na.rm = T),
          C_uneat = sum(sens[['C_uneat_stat']][,2], na.rm = T),
          food_prov = sum(sens[['food_prov_stat']][,2], na.rm = T),
          rel_feeding = mean(sens[['rel_feeding_stat']][,2], na.rm = T),
          ing_act = sum(sens[['ing_act_stat']][,2], na.rm = T),
          anab = mean(sens[['anab_stat']][,2], na.rm = T),
          catab = mean(sens[['catab_stat']][,2], na.rm = T),
          O2 = sum(sens[['O2_stat']][,2], na.rm = T),
          NH4 = sum(sens[['NH4_stat']][,2], na.rm = T),
          adj_param = as.factor(tar_param_names_pop),
          factor = fact,
          farm_ID = tar_farm_IDs
        )
      }) %>%
        pivot_longer(
          cols = c(weight, dw, total_excr, P_excr, L_excr, C_excr, total_uneat, P_uneat, L_uneat, C_uneat,
                   food_prov, rel_feeding, ing_act, anab, catab, O2, NH4),
          names_to = "measure", names_transform = list(measure = as.factor), values_to = "value"
        ) %>%
        group_by(adj_param, factor, farm_ID, measure) %>%
        pivot_wider(names_from = "factor", names_prefix = "fact_", values_from = "value") %>%
        mutate(sensitivity = (fact_1.1 - fact_0.9) / (0.2 * fact_1)) %>%
        group_by(adj_param, measure, farm_ID) %>%
        reframe(sensitivity = meanna(sensitivity))
    },
    pattern = cross(tar_param_names_pop, map(tar_farm_ts, tar_farm_IDs))
  ),
  
  tar_target(
    tar_sens_results_pop,
    command = {
      tar_sens_run_pop %>%
        group_by(adj_param, measure) %>%
        reframe(mean_sens = meanna(sensitivity),
                sd_sens = sdna(sensitivity))
    },
    pattern = tar_param_names_pop
  )
)
