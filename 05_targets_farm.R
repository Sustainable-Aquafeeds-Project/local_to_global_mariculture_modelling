suppressMessages(suppressWarnings(library(targets)))
suppressMessages(suppressWarnings(library(crew)))
suppressMessages(suppressWarnings(library(mirai)))
suppressMessages(suppressWarnings(library(arrow)))
suppressMessages(suppressWarnings(library(sf)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(terra))) # do not remove
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(conflicted)))
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "future", "furrr", "ggplot2", "matrixStats", "tibble"), 
  format = "qs", 
  controller = crew_controller_local(workers = 3),
  workspace_on_error = TRUE
)

tar_source("00_model_functions.R")

# Paths -----------------------------------------------------------------------------------------------------------
gendata_path <- file.path("data", "_general_data")
this_species <- "atlantic_salmon"
species_path <- file.path("data", this_species)

# Global objects --------------------------------------------------------------------------------------------------
farm_harvest_size <- tar_read(farm_harvest_size, store = "04_targets_individual") %>% 
  select(c(farm_ID, weight))

# Begin targets pipeline ------------------------------------------------------------------------------------------
list(
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
      slice_sample(n = 1125) %>%
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
  tar_target(farm_times, 
             c(farm_coords$t_start[farm_coords$farm_ID == farm_IDs], farm_coords$t_end[farm_coords$farm_ID == farm_IDs], dt = 1),
             pattern = farm_IDs),
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
      slice_head(n = 1),
    pattern = farm_IDs
  ),
  tar_target(
    N_population, 
    generate_pop(
      harvest_n = (farm_static_data$tonnes_per_farm*10^6)/farm_harvest_size$weight[farm_harvest_size$farm_ID == farm_IDs],
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
      N_pop = N_population/2,
      nruns = pop_params['nruns']
    ),
    pattern = cross(map(feed_types, feed_params_proteins, feed_params_carbs, feed_params_lipids), 
                    map(farm_IDs, farm_temp, farm_times, N_population))
  ),

  tar_target(CommSize, farm_harvest_size$weight[farm_harvest_size$farm_ID == farm_IDs], pattern = farm_IDs),
  tar_target(
    days_to_CommSize,
    data.frame(
      farm_ID = farm_IDs,
      feed = feed_types,
      Ub = days_to_CS(weight_stat = (main_farm_growth[[1]][,1] + main_farm_growth[[1]][,2]), days = 1:455, CS = CommSize),
      Mean = days_to_CS(weight_stat = main_farm_growth[[1]][,1], days = 1:455, CS = CommSize),
      Lb = days_to_CS(weight_stat = (main_farm_growth[[1]][,1] - main_farm_growth[[1]][,2]), days = 1:455, CS = CommSize)
    ),
    pattern = map(main_farm_growth, cross(feed_types, map(farm_IDs, CommSize))), 
    deployment = "main"
  ),
  
  tar_target(
    cohorts_biomass,
    command = {
     df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[2]])
     df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
     df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
     df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
     df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>% 
       mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_dw,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[3]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_SGR,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[4]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_P_excr,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[6]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_L_excr,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[7]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_C_excr,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[8]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_P_uneat,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[9]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_L_uneat,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[10]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_C_uneat,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[11]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_ing_act,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[12]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_O2,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[15]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_NH4,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[16]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_food_prov,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[17]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_total_excr_mat,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[20]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
    deployment = "main"
  ),
  tar_target(
    cohorts_total_uneat_mat,
    command = {
      df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[21]])
      df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
      df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
      df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
      df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>%
        mutate(farm_ID = farm_IDs, feed = feed_types)
    },
    pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
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



