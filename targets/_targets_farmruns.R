suppressPackageStartupMessages(suppressWarnings({
  library(targets)
  library(crew)
  library(qs)
  library(tidyverse)
}))

tar_option_set(
  format = "qs", 
  controller = crew::crew_controller_local(
    workers = 10, 
    seconds_idle = 120,
    crashes_max = 25
    ),
  workspace_on_error = TRUE,
  garbage_collection = 25,
  error = "abridge"
)

tar_source(
  files = c(
    here::here("src/dirs.R"),
    here::here("src/functions.R"),
    here::here("src/model_functions.R")
  )
)

# Dirs
output_path <- here::here("outputs")
output_farm_data_path <- file.path(output_path, "farm_data")

# Globals
# For quickly running a small sample, change inds_per_farm to a smaller number and uncomment "slice_sample()" line in farm_IDs below
farm_chunk_size <- 12
inds_per_farm <- 2500

list(
# Get previously saved data ---------------------------------------------------------------------------------------
  tar_target(species_params_file, file.path(output_species_data_path, "species_params.qs"), format = "file"),
  tar_target(pop_params_file, file.path(output_species_data_path, "pop_params.qs"), format = "file"),
  tar_target(farm_temporal_data_file, file.path(input_farm_coords_path, "farm_temporal_data.qs"), format = "file"),
  tar_target(farm_static_data_file, file.path(input_farm_coords_path, "farm_static_data_w_stocking.qs"), format = "file"),
  tar_target(harvest_size_file, file.path(gendata_path, "harvest_sizes", "farm_harvest_sizes.qs"), format = "file"),

  tar_target(feed_params_file_MD, file.path(output_species_data_path, "feed_params_MD.qs"), format = "file"),
  tar_target(feed_params_file_PD, file.path(output_species_data_path, "feed_params_PD.qs"), format = "file"),
  tar_target(feed_params_file_AD, file.path(output_species_data_path, "feed_params_AD.qs"), format = "file"),

  # tar_target(farm_coords, qread(farm_coords_file)),
  tar_target(
    farm_IDs,
    command = {
      farm_static_data_file %>% 
        qread() %>% 
        distinct(farm_id) %>% 
        # slice_head(n = 960) %>% 
        pull(farm_id)
    }
  ),

  tar_target(
    farm_IDs_chunked,
    split(farm_IDs, ceiling(seq_along(farm_IDs)/farm_chunk_size))
  ),

  tar_target(feed_params_MD, qread(feed_params_file_MD)), 
  tar_target(feed_names_MD, names(feed_params_MD)),
  tar_target(feed_params_PD, qread(feed_params_file_PD)), 
  tar_target(feed_names_PD, names(feed_params_PD)),
  tar_target(feed_params_AD, qread(feed_params_file_AD)), 
  tar_target(feed_names_AD, names(feed_params_AD)),

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
        feed_params = feed_params_PD[["plant_dominant_biomar"]],
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

  tar_target(all_params, c(qread(species_params_file), qread(pop_params_file))),

# Prepare parameters and farm data --------------------------------------------------------------------------------
  tar_target(
    harvest_size,
    command = {
      qread(harvest_size_file)
    },
    deployment = "main"
  ),

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

# Run growth ------------------------------------------------------------------------------------------------------
  tar_target(
    farm_run_MD_chunked,
    command = {# error == T
      ts_data <- farm_ts_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      static_data <- farm_static_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      purrr::map2_dfr(ts_data, static_data, function(ts, static) {
        fg <- farm_growth(
          pop_params = all_params,
          species_params = all_params,
          water_temp = ts$temp_c,
          feed_params = feed_params_MD[[feed_names_MD]],
          times = c(t_start = static$t_start, t_end = static$t_end, dt = 1),
          N_pop = ts$Npop,
          nruns = inds_per_farm
        )
        fg %>% 
          purrr::map_dfr(function(mat) {
            mat %>% 
              as.data.frame() %>% 
              rename(t = V1, mean = V2, sd = V3) %>% 
              mutate(
                farm_ID = as.integer(static$farm_ID),
                t = as.integer(t)
              )
          }, .id = "measure") %>% 
          mutate(measure = as.factor(measure))
      }) %>% 
      mutate(feed = as.factor(feed_names_MD))
    },
    pattern = cross(map(farm_ts_data, farm_static_data), feed_names_MD)
  ),

  tar_target(
    farm_run_PD_chunked,
    command = {# error == T
      ts_data <- farm_ts_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      static_data <- farm_static_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      purrr::map2_dfr(ts_data, static_data, function(ts, static) {
        fg <- farm_growth(
          pop_params = all_params,
          species_params = all_params,
          water_temp = ts$temp_c,
          feed_params = feed_params_PD[[feed_names_PD]],
          times = c(t_start = static$t_start, t_end = static$t_end, dt = 1),
          N_pop = ts$Npop,
          nruns = inds_per_farm
        )
        fg %>% 
          purrr::map_dfr(function(mat) {
            mat %>% 
              as.data.frame() %>% 
              rename(t = V1, mean = V2, sd = V3) %>% 
              mutate(
                farm_ID = as.integer(static$farm_ID),
                t = as.integer(t)
              )
          }, .id = "measure") %>% 
          mutate(measure = as.factor(measure))
      }) %>% 
      mutate(feed = as.factor(feed_names_PD))
    },
    pattern = cross(map(farm_ts_data, farm_static_data), feed_names_PD)
  ),

  tar_target(
    farm_run_AD_chunked,
    command = {# error == T
      ts_data <- farm_ts_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      static_data <- farm_static_data %>% 
        group_by(farm_ID) %>% 
        group_split()

      purrr::map2_dfr(ts_data, static_data, function(ts, static) {
        fg <- farm_growth(
          pop_params = all_params,
          species_params = all_params,
          water_temp = ts$temp_c,
          feed_params = feed_params_AD[[feed_names_AD]],
          times = c(t_start = static$t_start, t_end = static$t_end, dt = 1),
          N_pop = ts$Npop,
          nruns = inds_per_farm
        )
        fg %>% 
          purrr::map_dfr(function(mat) {
            mat %>% 
              as.data.frame() %>% 
              rename(t = V1, mean = V2, sd = V3) %>% 
              mutate(
                farm_ID = as.integer(static$farm_ID),
                t = as.integer(t)
              )
          }, .id = "measure") %>% 
          mutate(measure = as.factor(measure))
      }) %>% 
      mutate(feed = as.factor(feed_names_AD))
    },
    pattern = cross(map(farm_ts_data, farm_static_data), feed_names_AD)
  ),

# Process results ------------------------------------------------------------------------------------------------------
  tar_target(
    farm_results,
    command = {# error == T
      rbind(
        farm_run_MD_chunked %>% 
          group_by(farm_ID) %>% 
          group_split() %>% 
          lapply(function(farm_df) {
            farm_df %>% 
              mutate(prod_t = as.integer(t-min(t)+1))

            # Get biomass only (is already cumulative)
            biomass <- farm_df %>% 
              filter(measure == "biomass_stat" & t == max(t)) %>% 
              select(-t)

            # These stats are daily, need to accumulate
            full_df <- farm_df %>% 
              filter(measure %in% c("dw_stat", "food_prov_stat", "total_uneat_stat", "total_excr_stat", "L_excr_stat", "L_uneat_stat", "P_excr_stat", "P_uneat_stat", "C_excr_stat", "C_uneat_stat")) %>% 
              group_by(farm_ID, feed, measure) %>% 
              reframe(
                mean = sum(mean),
                sd = sqrt(sumna(sd^2))
              ) %>% 
              rbind(biomass) %>% 
              mutate(measure = droplevels(measure))
            
            # Calculate everything per biomass as well
            full_df %>% 
              left_join(
                biomass %>% 
                  rename(
                    mean_biomass = mean,
                    sd_biomass = sd
                  ) %>% select(-measure),
                  by = c("farm_ID", "feed")
              ) %>% 
              mutate(
                mean_per_biom = mean/mean_biomass,
                sd_per_biom = (mean/mean_biomass) * sqrt((sd_biomass/mean_biomass)^2 + (sd/mean)^2)
              )
            }) %>% 
          bind_rows(),
        farm_run_PD_chunked %>% 
          group_by(farm_ID) %>% 
          group_split() %>% 
          lapply(function(farm_df) {
            farm_df %>% 
              mutate(prod_t = as.integer(t-min(t)+1))

            # Get biomass only (is already cumulative)
            biomass <- farm_df %>% 
              filter(measure == "biomass_stat" & t == max(t)) %>% 
              select(-t)

            # These stats are daily, need to accumulate
            full_df <- farm_df %>% 
              filter(measure %in% c("dw_stat", "food_prov_stat", "total_uneat_stat", "total_excr_stat", "L_excr_stat", "L_uneat_stat", "P_excr_stat", "P_uneat_stat", "C_excr_stat", "C_uneat_stat")) %>% 
              group_by(farm_ID, feed, measure) %>% 
              reframe(
                mean = sum(mean),
                sd = sqrt(sumna(sd^2))
              ) %>% 
              rbind(biomass) %>% 
              mutate(measure = droplevels(measure))
            
            # Calculate everything per biomass as well
            full_df %>% 
              left_join(
                biomass %>% 
                  rename(
                    mean_biomass = mean,
                    sd_biomass = sd
                  ) %>% select(-measure),
                  by = c("farm_ID", "feed")
              ) %>% 
              mutate(
                mean_per_biom = mean/mean_biomass,
                sd_per_biom = (mean/mean_biomass) * sqrt((sd_biomass/mean_biomass)^2 + (sd/mean)^2)
              )
            }) %>% 
          bind_rows(),
        farm_run_AD_chunked %>% 
          group_by(farm_ID) %>% 
          group_split() %>% 
          lapply(function(farm_df) {
            farm_df %>% 
              mutate(prod_t = as.integer(t-min(t)+1))

            # Get biomass only (is already cumulative)
            biomass <- farm_df %>% 
              filter(measure == "biomass_stat" & t == max(t)) %>% 
              select(-t)

            # These stats are daily, need to accumulate
            full_df <- farm_df %>% 
              filter(measure %in% c("dw_stat", "food_prov_stat", "total_uneat_stat", "total_excr_stat", "L_excr_stat", "L_uneat_stat", "P_excr_stat", "P_uneat_stat", "C_excr_stat", "C_uneat_stat")) %>% 
              group_by(farm_ID, feed, measure) %>% 
              reframe(
                mean = sum(mean),
                sd = sqrt(sumna(sd^2))
              ) %>% 
              rbind(biomass) %>% 
              mutate(measure = droplevels(measure))
            
            # Calculate everything per biomass as well
            full_df %>% 
              left_join(
                biomass %>% 
                  rename(
                    mean_biomass = mean,
                    sd_biomass = sd
                  ) %>% select(-measure),
                  by = c("farm_ID", "feed")
              ) %>% 
              mutate(
                mean_per_biom = mean/mean_biomass,
                sd_per_biom = (mean/mean_biomass) * sqrt((sd_biomass/mean_biomass)^2 + (sd/mean)^2)
              )
            }) %>% 
          bind_rows()
      )
    },
    pattern = map(farm_run_MD_chunked, farm_run_PD_chunked, farm_run_AD_chunked),
    deployment = "main"
  ),

  tar_target(
    total_vals,
    command = {
      # Some initial cleanup
      total_vals_0 <- farm_results %>% 
        mutate(
          feed_group = factor(
            str_remove(feed, "_min|_max"), 
            levels = c("marine_dominant_biomar", "plant_dominant_biomar", "animal_dominant_biomar"),
            labels = c("MD", "PD", "AD")
            ),
          feed_type = case_when(
            str_detect(feed, "_min") ~ "min",
            str_detect(feed, "_max") ~ "max",
            T ~ "mean"
            ),
          feed_type = factor(
            feed_type, 
            levels = c("min", "mean", "max"), 
            labels = c("Minimum", "Mean", "Maximum")
          )
        ) %>% 
        select(-feed)

      # Sum some measures - mean_per_biom
      total_vals_1 <- total_vals_0 %>% 
        filter(measure %in% c("food_prov_stat", "total_excr_stat", "total_uneat_stat", "C_excr_stat", "C_uneat_stat","L_excr_stat", "L_uneat_stat","P_excr_stat", "P_uneat_stat")) %>% 
        mutate(measure = droplevels(measure)) %>% 
        select(farm_ID, measure, mean_per_biom, sd_per_biom, feed_group, feed_type) %>% 
        pivot_wider(
          names_from = measure, names_sep = ".",
          values_from = c(mean_per_biom, sd_per_biom)
        ) %>% 
        mutate(
          mean_per_biom.total_waste_stat = mean_per_biom.total_excr_stat + mean_per_biom.total_uneat_stat,
          sd_per_biom.total_waste_stat = mean_per_biom.total_waste_stat * (sd_per_biom.total_excr_stat/mean_per_biom.total_excr_stat + sd_per_biom.total_uneat_stat/mean_per_biom.total_uneat_stat),
          mean_per_biom.total_P_lost_stat = mean_per_biom.P_excr_stat + mean_per_biom.P_uneat_stat,
          sd_per_biom.total_P_lost_stat = mean_per_biom.total_P_lost_stat * (sd_per_biom.P_excr_stat/mean_per_biom.P_excr_stat + sd_per_biom.P_uneat_stat/mean_per_biom.P_uneat_stat),
          mean_per_biom.total_C_lost_stat = mean_per_biom.C_excr_stat + mean_per_biom.C_uneat_stat,
          sd_per_biom.total_C_lost_stat = mean_per_biom.total_C_lost_stat * (sd_per_biom.C_excr_stat/mean_per_biom.C_excr_stat + sd_per_biom.C_uneat_stat/mean_per_biom.C_uneat_stat),
          mean_per_biom.total_L_lost_stat = mean_per_biom.L_excr_stat + mean_per_biom.L_uneat_stat,
          sd_per_biom.total_L_lost_stat = mean_per_biom.total_L_lost_stat * (sd_per_biom.L_excr_stat/mean_per_biom.L_excr_stat + sd_per_biom.L_uneat_stat/mean_per_biom.L_uneat_stat),
        ) %>% 
        pivot_longer(
          cols = contains("stat"),
          names_to = c(".value", "measure"),
          names_transform = list(measure = as.factor),
          names_sep = "\\.",
        )

      # Sum some measures - mean
      total_vals_2 <- total_vals_0 %>% 
        filter(measure %in% c("food_prov_stat", "total_excr_stat", "total_uneat_stat", "C_excr_stat", "C_uneat_stat", "L_excr_stat", "L_uneat_stat", "P_excr_stat", "P_uneat_stat")) %>% 
        mutate(measure = droplevels(measure)) %>% 
        select(farm_ID, measure, mean, sd, feed_group, feed_type) %>% 
        pivot_wider(
          names_from = measure, 
          values_from = c(mean, sd),
          names_sep = "."
        ) %>% 
        mutate(
          mean.total_waste_stat = mean.total_excr_stat + mean.total_uneat_stat,
          sd.total_waste_stat = mean.total_waste_stat * (sd.total_excr_stat/mean.total_excr_stat + sd.total_uneat_stat/mean.total_uneat_stat),
          mean.total_P_lost_stat = mean.P_excr_stat + mean.P_uneat_stat,
          sd.total_P_lost_stat = mean.total_P_lost_stat * (sd.P_excr_stat/mean.P_excr_stat + sd.P_uneat_stat/mean.P_uneat_stat),
          mean.total_C_lost_stat = mean.C_excr_stat + mean.C_uneat_stat,
          sd.total_C_lost_stat = mean.total_C_lost_stat * (sd.C_excr_stat/mean.C_excr_stat + sd.C_uneat_stat/mean.C_uneat_stat),
          mean.total_L_lost_stat = mean.L_excr_stat + mean.L_uneat_stat,
          sd.total_L_lost_stat = mean.total_L_lost_stat * (sd.L_excr_stat/mean.L_excr_stat + sd.L_uneat_stat/mean.L_uneat_stat)
        ) %>% 
        pivot_longer(
          cols = contains("stat"),
          names_to = c(".value", "measure"),
          names_transform = list(measure = as.factor),
          names_sep = "\\.",
        )

      # Rejoin mean and mean_per_biom
      left_join(total_vals_1, total_vals_2, by = c("farm_ID", "feed_group", "feed_type", "measure"))
    },
    pattern = farm_results,
    deployment = "main"
  ),

  tar_target(
    total_props,
    command = {
      total_vals %>% 
        pivot_wider(
          names_from = measure, names_sep = ".",
          values_from = c(mean, mean_per_biom)
        ) %>% 
        mutate(
          prop_excr_stat = mean.total_excr_stat/mean.total_waste_stat, # The proportion of total waste (feed+faeces) that comes from excreted faeces
          prop_uneat_stat = mean.total_uneat_stat/mean.food_prov_stat, # The proportion of feed provided that is lost to the environment
          prop_P_stat = mean.P_excr_stat/mean.total_P_lost_stat, # The proportion of total protein lost (feed+faeces) that is lost through excretion
          prop_C_stat = mean.C_excr_stat/mean.total_C_lost_stat, # The proportion of total carbohydrate lost (feed+faeces) that is lost through excretion
          prop_L_stat = mean.L_excr_stat/mean.total_L_lost_stat, # The proportion of total lipid lost (feed+faeces) that is lost through excretion
          prop_comp_excr_P_stat = mean.P_excr_stat/mean.total_excr_stat, # The composition of faeces, ie how much of faeces is protein
          prop_comp_excr_C_stat = mean.C_excr_stat/mean.total_excr_stat, # The composition of faeces, ie how much of faeces is carbohydrate
          prop_comp_excr_L_stat = mean.L_excr_stat/mean.total_excr_stat, # The composition of faeces, ie how much of faeces is lipid
        ) %>% 
        select(farm_ID, feed_group, feed_type, contains("prop")) %>% 
        pivot_longer(
          cols = contains("prop"), 
          names_to = c("measure"), 
          names_transform = as.factor,
          values_to = "mean_prop", 
        )
    },
    pattern = total_vals,
    deployment = "main"
  )

  # tar_target(
  #   cohort_results,
  #   command = {
  #     farm_results <- split(farm_run_chunked, farm_run_chunked$farm_ID) %>% 
  #       lapply(function(farm_df) {
  #         # Define the three cohorts
  #         cohort_1 <- farm_df %>% 
  #           mutate(prod_t = as.integer(t-min(t)+1)) %>% 
  #           filter(
  #             prod_t != max(prod_t) & 
  #               !measure %in% c("anab_stat", "catab_stat","E_assim_stat", "E_somat_stat", "food_enc_stat", "ing_act_stat", "ing_pot_stat", "metab_stat", "NH4_stat", "O2_stat", "rel_feeding_stat", "T_response_stat", "water_temp_stat", "weight_stat")
  #           ) %>% 
  #           mutate(
  #             measure = droplevels(measure),
  #             cohort = 1
  #           )
  #         cohort_2 <- farm_to_cohort(cohort_1, time_offset = 365) %>% mutate(cohort = 2)
  #         cohort_3 <- farm_to_cohort(cohort_1, time_offset = 730) %>% mutate(cohort = 3)

  #         # Start and end times for each
  #         starts <- rbind(
  #           cohort_2 %>% distinct(cohort, t, prod_t) %>% slice_min(t, n = 1),
  #           cohort_3 %>% distinct(cohort, t, prod_t) %>% slice_min(t, n = 1)
  #         )
  #         ends <- rbind(
  #           cohort_1 %>% distinct(cohort, t, prod_t) %>% slice_max(t, n = 1),
  #           cohort_2 %>% distinct(cohort, t, prod_t) %>% slice_max(t, n = 1)
  #         )

  #         # Calculate stats based on daily values across overlapping cohorts
  #         full_df <- rbind(cohort_1, cohort_2, cohort_3) %>% 
  #           filter(t %in% starts$t[starts$cohort == 2]:ends$t[ends$cohort == 2]) %>% 
  #           group_by(farm_ID, feed, measure, t) %>%
  #           reframe(
  #             mean = sumna(mean),
  #             sd = sqrt(sumna(sd^2))) %>% 
  #           mutate(
  #             cohorts = as.factor(case_when(
  #               t <= ends$t[ends$cohort == 1] ~ "1 & 2",
  #               t >= starts$t[starts$cohort == 3] ~ "2 & 3",
  #               T ~ "2"
  #             )),
  #             prod_t = as.integer(t-min(t)+1)
  #           )
          
  #         # Combine with biomass to get values per biomass produced
  #         full_df %>% 
  #           left_join(
  #             full_df %>% filter(measure == "biomass_stat") %>% select(-measure),
  #             by = c("farm_ID", "cohorts", "feed", "t", "prod_t")
  #           ) %>% 
  #           rename(
  #             mean = mean.x, 
  #             sd = sd.x,
  #             mean_biomass = mean.y, 
  #             sd_biomass = sd.y
  #           ) %>% 
  #           mutate(
  #             mean_per_biom = mean/mean_biomass,
  #             sd_per_biom = (mean/mean_biomass) * sqrt((sd_biomass/mean_biomass)^2 + (sd/mean)^2)
  #           )
  #       }) %>% 
  #       bind_rows()
  #   },
  #   pattern = farm_run_chunked,
  #   deployment = "main"
  # )
)
