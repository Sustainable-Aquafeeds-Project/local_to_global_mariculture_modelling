library(targets)
library(crew)

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "tibble", "qs"), 
  format = "qs", 
  controller = crew_controller_local(workers = 10),
  workspace_on_error = TRUE
)

tar_source(files = list.files("src", pattern = "\\.R$", full.names = TRUE))

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
      # c(1:39, 2704:2721) # Australia and US only

      qs::qread(farm_ts_data_file) %>% 
        distinct(farm_ID) %>% 
        pull(farm_ID)
    }
  ),
  tar_target(
    farm_ts_data, 
    qs::qread(farm_ts_data_file) %>% 
      merge(farm_coords, by = "farm_ID") %>%
      filter(farm_ID == farm_IDs) %>%
      mutate(keep = case_when(day >= t_start & day <= t_end ~ T, T ~ F)) %>% 
      dplyr::filter(keep) %>% 
      dplyr::select(farm_ID, day, temp_c),
    pattern = farm_IDs
  ),
  tar_target(
    feed_names,
    command = {
      qs::qread(feed_params_file) %>% names()
      # c("plant_dominant", "marine_dominant")
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
      farm_temp <- farm_ts_data %>%
        filter(farm_ID == fi)
      fg <- fish_growth(
        pop_params = all_params,
        species_params = all_params,
        water_temp = farm_temp$temp_c,
        feed_params = reference_feed,
        times = c(t_start = min(farm_temp$day), t_end = max(farm_temp$day), dt = 1),
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
    harvest_size,
    command = {
      fg <- fish_growth(
        pop_params = all_params,
        species_params = all_params,
        water_temp = farm_ts_data$temp_c,
        feed_params = reference_feed,
        times = c(t_start = min(farm_ts_data$day), t_end = max(farm_ts_data$day), dt = 1),
        init_weight = all_params["meanW"],
        ingmax = all_params["meanImax"]
      ) 
      fg %>%
        as.data.frame() %>%
        filter(!is.na(weight)) %>% 
        slice_tail(n = 1) %>%
        select(weight) %>%
        mutate(
          farm_ID = farm_IDs,
          weight = weight %>% 
            units::set_units("g") %>% 
            units::drop_units()
        )
    },
    pattern = map(farm_IDs, farm_ts_data)
  ),

# Prepare parameters ----------------------------------------------------------------------------------------------
  tar_target(
    stat_names,
    c(
      "weight_stat", "dw_stat", "water_temp_stat", "T_response_stat", "P_excr_stat", "L_excr_stat", "C_excr_stat", 
      "P_uneat_stat", "L_uneat_stat", "C_uneat_stat", "food_prov_stat", "food_enc_stat", "rel_feeding_stat", 
      "ing_pot_stat", "ing_act_stat", "E_assim_stat", "E_somat_stat", "anab_stat","catab_stat", "O2_stat", 
      "NH4_stat", "total_excr_stat", "total_uneat_stat", "metab_stat", "biomass_stat"
    ),
  ),

  tar_target(
    farm_ID_data, 
    command = {
      farm_production <- farm_production_file %>% 
        qs::qread() %>% 
        distinct(farm_id, tonnes_per_farm) %>% 
        mutate(
          tonnes_per_farm = tonnes_per_farm %>% 
            units::set_units("t") %>% 
            units::set_units("g") %>% 
            units::drop_units()
        )
      c(
        't_start' = min(farm_ts_data$day),
        't_end' = max(farm_ts_data$day),
        'dt' = 1,
        'nruns' = 5000,
        'prod' = farm_production$tonnes_per_farm[farm_production$farm_id == farm_IDs]
      )
    },
    pattern = map(farm_IDs, farm_ts_data)
  ),
  
  tar_target(
    N_pop, 
    command = {
      generate_pop(
        harvest_n = farm_ID_data['prod']/harvest_size$weight,
        mort = all_params['mortmyt'],
        times = farm_ID_data
      )
    }, 
    pattern = map(farm_IDs, harvest_size, farm_ID_data)
  ),
  
# Run growth ------------------------------------------------------------------------------------------------------
  tar_target(
    farm_run,
    command = {
      fg <- farm_growth(
        pop_params = all_params,
        species_params = all_params,
        water_temp = farm_ts_data$temp_c,
        feed_params = feed_params,
        times = farm_ID_data,
        N_pop = N_pop,
        nruns = farm_ID_data['nruns']
      ) 
      fg %>% 
        lapply(., function(mat) {
          mat %>% 
            as.data.frame() %>% 
            rename(t = V1, mean = V2, sd = V3) %>% 
            mutate(
              farm_ID = farm_IDs,
              feed = as.factor(feed_names)
            )
        })
    },
    pattern = cross(map(farm_IDs, N_pop, farm_ID_data, farm_ts_data), map(feed_params, feed_names))
  ),

  tar_target(
    farm_results,
    command = {
      farm_run[[stat_names]] %>% 
        mutate(sd_mean = sd/mean) %>%
        group_by(farm_ID, feed, t) %>%
        reframe(
          mean = sumna(mean),
          sd_mean = sumna(sd_mean),
          sd = sd_mean*mean,
          measure = as.factor(stat_names)
        ) %>% 
      # Fish are given no food on harvest day
      group_by(farm_ID) %>% 
      mutate(last_t = max(t)) %>% 
      ungroup() %>% 
      filter(t != last_t) %>% 
      dplyr::select(-last_t)
    },
    pattern = cross(farm_run, stat_names)
  ),

  tar_target(
    cohort_results,
    command = {
      cohort_1 <- farm_run[[stat_names]] %>% 
        mutate(sd_mean = sd/mean) %>% 
      # Fish are given no food on harvest day
      group_by(farm_ID) %>% 
      mutate(last_t = max(t)) %>% 
      ungroup() %>% 
      filter(t != last_t) %>% 
      dplyr::select(-last_t)
      
      cohort_2 <- farm_to_cohort(cohort_1, time_offset = 365)
      cohort_3 <- farm_to_cohort(cohort_1, time_offset = 730)
      cohort_1 <- farm_to_cohort(cohort_1, time_offset = 0)
      lims <- unique(cohort_2$t)

      rbind(cohort_1, cohort_2, cohort_3) %>% 
        filter(t %in% lims) %>% 
        group_by(farm_ID, feed, t) %>%
        reframe(
          mean = sumna(mean),
          sd_mean = sumna(sd_mean),
          sd = sd_mean*mean,
          measure = as.factor(stat_names)
        )
    },
    pattern = cross(farm_run, stat_names)
  )
)
