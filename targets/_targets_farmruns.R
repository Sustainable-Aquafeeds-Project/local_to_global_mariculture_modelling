suppressPackageStartupMessages(suppressWarnings({
  library(targets)
  library(crew)
  library(magrittr)
  library(qs)
  library(tidyverse)
  library(arrow)
}))


tar_option_set(
  # packages = c("arrow"), 
  format = "qs", 
  controller = crew::crew_controller_local(workers = 5, seconds_idle = 300),
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
        # slice_head(n = 24) %>% 
        pull(farm_ID)
    }
  ),

  tar_target(
    farm_IDs_chunked,
    split(farm_IDs, ceiling(seq_along(farm_IDs)/10))
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
    command = {names(qs::qread(feed_params_file))}
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

# Prepare parameters and farm data --------------------------------------------------------------------------------
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

# Run growth ------------------------------------------------------------------------------------------------------
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
          nruns = 500
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
    pattern = cross(map(farm_ts_data_chunked, farm_static_data_chunked), map(feed_params, feed_names))
  ),

# Process results ------------------------------------------------------------------------------------------------------
  tar_target(
    farm_full_results,
    command = {
      farm_results <- farm_run_chunked %>% 
        mutate(prod_t = t-min(t)+1)

      total_ends <- rbind(
        farm_results %>% 
          filter(measure == "biomass_stat" & t == max(t)) %>% 
          select(-t, -prod_t), 
        farm_results %>% 
          filter(measure %in% c("food_prov_stat", "total_uneat_stat", "total_excr_stat", "L_excr_stat", "L_uneat_stat", "P_excr_stat", "P_uneat_stat", "C_excr_stat", "C_uneat_stat")) %>% 
          group_by(farm_ID, feed, measure) %>% 
          reframe(
            mean = sum(mean),
            sd = sqrt(sumna(sd^2))
          )
      )
      total_ends <- merge(
        total_ends %>% filter(measure != "biomass_stat"),
        total_ends %>% filter(measure == "biomass_stat") %>% select(-measure),
        by = c("farm_ID", "feed")
        ) %>% 
        rename(
          mean = mean.x, sd = sd.x,
          mean_biomass = mean.y, sd_biomass = sd.y
        ) %>% 
        mutate(
          mean_per_biom = mean/mean_biomass,
          sd_per_biom = (mean/mean_biomass) * sqrt((sd_biomass/mean_biomass)^2 + (sd/mean)^2)
        )

      # Calculate stats based on daily values across cohorts
      cohort_results_daily <- farm_results %>% 
        filter(prod_t != max(prod_t) & !measure %in% c("anab_stat", "catab_stat", "dw_stat", "E_assim_stat", "E_somat_stat", "food_enc_stat", "ing_act_stat", "ing_pot_stat", "metab_stat", "NH4_stat", "O2_stat", "rel_feeding_stat", "T_response_stat", "water_temp_stat", "weight_stat"))
      cohort_results_daily <- split(cohort_results_daily, cohort_results_daily$farm_ID)
      cohort_results_daily <- purrr::map(cohort_results_daily, function(crd) {
        cohort_1 <- farm_to_cohort(crd, time_offset = 0)
        cohort_2 <- farm_to_cohort(crd, time_offset = 365)
        cohort_3 <- farm_to_cohort(crd, time_offset = 730)
        lims <- min(cohort_2$t):(min(cohort_2$t)+365*2) # Two years after start of cohort 2
        rbind(cohort_1, cohort_2, cohort_3) %>% 
          filter(t %in% lims) %>% 
          group_by(farm_ID, feed, measure, t) %>%
          reframe(
            mean = sumna(mean),
            sd = sqrt(sumna(sd^2))
          ) %>% 
          mutate(prod_t = t-min(lims)+1)
      }) %>% 
        bind_rows()

      cohort_results_daily <- merge(
        cohort_results_daily %>% filter(measure != "biomass_stat"),
        cohort_results_daily %>% filter(measure == "biomass_stat") %>% select(-measure),
        by = c("farm_ID", "feed", "t", "prod_t")
        ) %>% 
        rename(
          mean = mean.x, sd = sd.x,
          mean_biomass = mean.y, sd_biomass = sd.y
        ) %>% 
        mutate(
          mean_per_biom = mean/mean_biomass,
          sd_per_biom = (mean/mean_biomass) * sqrt((sd_biomass/mean_biomass)^2 + (sd/mean)^2)
        )
      
      list(
        total_ends = total_ends,
        cohort_results_daily = cohort_results_daily
      )
    },
    pattern = farm_run_chunked
  )

)
