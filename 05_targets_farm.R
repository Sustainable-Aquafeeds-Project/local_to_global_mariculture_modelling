suppressMessages(suppressWarnings(suppressPackageStartupMessages({
  library(targets)
  library(crew)
  library(mirai)
  library(arrow)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(terra) # do not remove
  library(magrittr)
  library(conflicted)
  conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)
})))

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "future", "furrr", "ggplot2", "matrixStats", "tibble"), 
  format = "qs", 
  controller = crew_controller_local(workers = parallelly::availableCores()-2),
  workspace_on_error = TRUE
)

tar_source("00_model_functions.R")

# Paths -----------------------------------------------------------------------------------------------------------
gendata_path <- file.path("data", "_general_data")
this_species <- "atlantic_salmon"
species_path <- file.path("data", this_species)

# Global objects --------------------------------------------------------------------------------------------------
farm_harvest_size <- tar_read(farm_harvest_size, store = "04_targets_individual") %>% 
  select(c(farm_ID, weight)) %>% 
  mutate(weight = units::set_units(weight, "g"))

# Begin targets pipeline ------------------------------------------------------------------------------------------
list(
  ## Farm data -----------------------------------------------------------------------------------------------------
  tar_target(farm_data_file, file.path(gendata_path, "SST", "farm_SST_extracted.parquet"), format = "file"),
  tar_target(farm_coord_file, file.path(gendata_path, "farm_locations", "farm_coords.parquet"), format = "file"),
  tar_target(farms_to_omit_file, sprintf(file.path(gendata_path, "farm_locations", "%s_farms_to_omit.qs"), this_species), format = "file"),
  tar_target(farms_to_omit, qs::qread(farms_to_omit_file)),
  tar_target(
    farm_IDs, 
    farm_data_file %>% 
      read_parquet() %>% 
      filter(!farm_id %in% farms_to_omit) %>% 
      distinct(farm_id) %>% 
      # slice_sample(n = 350) %>%
      as.vector() %>% unlist() %>% unname()
  ),
  tar_target(
    farm_ts_data, 
    read_parquet(farm_data_file) %>% 
      select(c(farm_id, day, temp_c)) %>% 
      mutate(day = str_split_i(day, "day_", 2) %>% as.numeric()) %>% 
      filter(farm_id == farm_IDs),
    pattern = farm_IDs
  ),
  
  # starts 1 May in northern hemisphere, 1 October in southern hemisphere
  tar_target(times_N, c("t_start" = 121, "t_end" = 121+547, "dt" = 1)),
  tar_target(times_S, c("t_start" = 274, "t_end" = 274+547, "dt" = 1)),
  tar_target(
    farm_coords, 
    read_parquet(farm_coord_file) %>% 
      mutate(t_start = case_when(lat > 0 ~ times_N['t_start'], TRUE ~ times_S['t_start']),
             t_end = case_when(lat > 0 ~ times_N['t_end'], TRUE ~ times_S['t_end']))
  ),
  tar_target(
    farm_times, 
    c(farm_coords$t_start[farm_coords$farm_ID == farm_IDs], farm_coords$t_end[farm_coords$farm_ID == farm_IDs], dt = 1),
    pattern = farm_IDs
  ),
  tar_target(farm_temp, farm_ts_data$temp_c[farm_times['t_start']:farm_times['t_end']], pattern = map(farm_ts_data, farm_times)),
  
  ## Species parameters -------------------------------------------------------------------------------------------
  tar_target(species_param_file, file.path(species_path, "params", "Species.xlsx"), format = "file"),
  tar_target(species_params, command = {
    df <- readxl::read_excel(species_param_file, sheet = "Atlantic salmon")
    vc <- df$Value
    names(vc) <- df$Quantity
    vc[!is.na(vc)]
  }),
  
  ## Population parameters ----------------------------------------------------------------------------------------
  tar_target(pop_params_file, file.path(species_path, "params", "Population.xlsx"), format = "file"),
  tar_target(pop_params, command = {
    df <- readxl::read_excel(pop_params_file)
    vc <- df$Value
    names(vc) <- df$Quantity
    vc[!is.na(vc)]
  }),
  tar_target(
    farm_static_data, 
    read_parquet(farm_data_file) %>% 
      select(-c(day, temp_c)) %>% 
      filter(farm_id == farm_IDs) %>% 
      slice_head(n = 1) %>% 
      mutate(tonnes_per_farm = tonnes_per_farm %>% units::set_units("t") %>% units::set_units("g")),
    pattern = farm_IDs
  ),
  tar_target(
    N_population, 
    generate_pop(
      harvest_n = 1.117811 * units::drop_units(farm_static_data$tonnes_per_farm/farm_harvest_size$weight[farm_harvest_size$farm_ID == farm_IDs]),
      mort = pop_params['mortmyt'],
      times = farm_times
    ), 
    pattern = map(farm_IDs, farm_static_data, farm_times)
  ),
  
  # Feed parameters -----------------------------------------------------------------------------------------------
  tar_target(feed_types, c("reference", "past", "future")),
  tar_target(feed_profile_file, file.path(gendata_path, "diets", "feed_profiles_and_composition_Atlantic_salmon.xlsx"), format = "file"),
  tar_target(
    feed_params_proteins, 
    readxl::read_excel(feed_profile_file, sheet = feed_types) %>% 
      select(ingredient, proportion, contains("protein")) %>% 
      select(-contains("feed")) %>% 
      rename(macro = ing_protein, digest = ing_protein_digestibility), 
    pattern = feed_types,
    iteration = "list"
  ),
  tar_target(
    feed_params_carbs, 
    readxl::read_excel(feed_profile_file, sheet = feed_types) %>% 
      select(ingredient, proportion, contains("carb")) %>% 
      select(-contains("feed")) %>% 
      rename(macro = ing_carb, digest = ing_carb_digestibility), 
    pattern = feed_types,
    iteration = "list"
  ),
  tar_target(
    feed_params_lipids, 
    readxl::read_excel(feed_profile_file, sheet = feed_types) %>% 
      select(ingredient, proportion, contains("lipid")) %>% 
      select(-contains("feed")) %>% 
      rename(macro = ing_lipid, digest = ing_lipid_digestibility), 
    pattern = feed_types,
    iteration = "list"
  ),
  
  # Main farm growth ----------------------------------------------------------------------------------------------
  tar_target(
    main_farm_growth,
    farm_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = farm_temp,
      feed_params = list(
        Proteins = feed_params_proteins,
        Carbohydrates = feed_params_carbs,
        Lipids = feed_params_lipids
      ),
      times = farm_times,
      N_pop = N_population,
      nruns = pop_params['nruns']
    ),
    pattern = cross(map(feed_types, feed_params_proteins, feed_params_carbs, feed_params_lipids), 
                    map(farm_IDs, farm_temp, farm_times, N_population))
  ),
  
  # [1] "main_farm_growth_weight_stat"
  # [2] "main_farm_growth_biomass_stat"
  # [3] "main_farm_growth_dw_stat"
  # [4] "main_farm_growth_SGR_stat"
  # [5] "main_farm_growth_E_somat_stat"
  # [6] "main_farm_growth_P_excr_stat"
  # [7] "main_farm_growth_L_excr_stat"
  # [8] "main_farm_growth_C_excr_stat"
  # [9] "main_farm_growth_P_uneat_stat"
  # [10] "main_farm_growth_L_uneat_stat"
  # [11] "main_farm_growth_C_uneat_stat"
  # [12] "main_farm_growth_ing_act_stat"
  # [13] "main_farm_growth_anab_stat"
  # [14] "main_farm_growth_catab_stat"
  # [15] "main_farm_growth_NH4_stat"
  # [16] "main_farm_growth_O2_stat"
  # [17] "main_farm_growth_food_prov_stat"
  # [18] "main_farm_growth_rel_feeding_stat"
  # [19] "main_farm_growth_T_response_stat"
  
  ## Cohorts data -------------------------------------------------------------------------------------------------
  # All this does is take main_farm_growth and overlay it onto itself 3 times for 3 different cohorts
  tar_target(
    cohorts_data,
    command = {
      cohorts_data <- list()
      df0 <- matrix(farm_times["t_start"]:farm_times["t_end"], 548, 1)
      stat_names <- c("weight", "biomass", "dw", "SGR", "E_somat", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "ing_act", "anab", "catab", "NH4", "O2", "food_prov", "rel_feeding", "T_response")
      for (i in c(2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17)) {
        df1 <- df2 <- df3 <- cbind(df0, main_farm_growth[[i]], matrix(1, nrow(df0), 1), matrix(1:nrow(df0), nrow(df0), 1))
        colnames(df1) <- colnames(df2) <- colnames(df3) <- c("yday", "mean", "sd", "cohort", "prod_day")
        df2[,'cohort'] <- 2
        df2[,'yday'] <- df2[,'yday']+365
        df3[,'cohort'] <- 3
        df3[,'yday'] <- df3[,'yday']+365+365
        cohorts_data[[i]] <- rbind(df1, df2, df3) %>% 
          as.data.frame() %>% tibble::remove_rownames() %>% 
          mutate(farm_ID = farm_IDs, feed = feed_types, stat = stat_names[i])
      }
      cohorts_data[!sapply(cohorts_data, is.null)]
    },
    pattern = map(main_farm_growth, cross(feed_types, map(farm_IDs, farm_times))), 
    deployment = "main"
  )
)

# library(units)
# df <- tar_read(cohorts_biomass, branches = 1)
# df$mean <- df$mean %>% set_units("g") %>% set_units("kg") %>% drop_units()
# whole_cohort <- df %>% filter(days %in% df$days[df$cohort == 2])
# whole_year <- df %>% mutate(days = days-365) %>% filter(days > 0 & days <= 365)
# 
# library(ggplot2)
# ggplot(whole_cohort, aes(x = days, y = mean, colour = as.factor(cohort))) +
#   geom_line(linewidth = 0.75) +
#   theme_classic()
# ggplot(whole_year, aes(x = days, y = mean, colour = as.factor(cohort))) +
#   geom_line(linewidth = 0.75) +
#   theme_classic()



