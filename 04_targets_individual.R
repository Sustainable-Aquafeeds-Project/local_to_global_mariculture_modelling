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
  library(units)
  conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)
})))

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "future", "furrr", "ggplot2", "matrixStats", "tibble"), 
  format = "qs", 
  controller = crew_controller_local(workers = parallelly::availableCores()-2, seconds_idle = 10),
  workspace_on_error = TRUE
)

tar_source("00_model_functions.R")

# Paths -----------------------------------------------------------------------------------------------------------
this_species <- "atlantic_salmon"
this_path <- file.path("data", this_species)
google_path <- "C:/Users/treimer/OneDrive - University of Tasmania/Sustainable Aquafeeds/Data from Google Drive"

# Begin targets pipeline ------------------------------------------------------------------------------------------
list(
  tar_target(farm_data_file, "data/_general_data/SST/farm_SST_extracted.parquet", format = "file"),
  tar_target(farm_coord_file, "data/_general_data/farm_locations/farm_coords.parquet", format = "file"),
  tar_target(farms_to_omit_file, sprintf("data/_general_data/farm_locations/%s_farms_to_omit.qs", this_species), format = "file"),
  tar_target(farms_to_omit, qs::qread(farms_to_omit_file)),
  tar_target(
    farm_IDs, 
    farm_data_file %>% 
      read_parquet() %>% 
      filter(!farm_id %in% farms_to_omit) %>% 
      distinct(farm_id) %>% 
      # slice_sample(n = 272) %>% 
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
  
  # Species parameters
  tar_target(species_param_file, file.path(this_path, "params", "Species.xlsx"), format = "file"),
  tar_target(species_params, command = {
    df <- readxl::read_excel(species_param_file, sheet = "Atlantic salmon")
    vc <- df$Value
    names(vc) <- df$Quantity
    vc[!is.na(vc)]
  }),
  
  # Population parameters
  tar_target(pop_params_file, file.path(this_path, "params", "Population.xlsx"), format = "file"),
  tar_target(pop_params, command = {
    df <- readxl::read_excel(pop_params_file)
    vc <- df$Value
    names(vc) <- df$Quantity
    vc[!is.na(vc)]
  }),
  
  # Feed parameters
  tar_target(feed_types, c("reference", "past", "future")),
  tar_target(feed_profile_file, "data/_general_data/diets/feed_profiles_and_composition_Atlantic_salmon.xlsx", format = "file"),
  tar_target(
    feed_params_protein, 
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
  
  ## Example individuals --------------------------------------------------------------------------------------------
  # tar_target(
  #   example_individual,
  #   command = {
  #     df <- fish_growth(
  #       pop_params = pop_params,
  #       species_params = species_params,
  #       water_temp = farm_temp,
  #       feed_params = list(
  #         Proteins = feed_params_protein,
  #         Carbohydrates = feed_params_carbs,
  #         Lipids = feed_params_lipids
  #       ),
  #       times = farm_times,
  #       init_weight = pop_params["meanW"],
  #       ingmax = pop_params["meanImax"]
  #     )
  #     df %>% 
  #       as.data.frame() %>% remove_rownames() %>% 
  #       mutate(SGR = 100 * (exp((log(weight)-log(weight[1]))/(days-days[1])) - 1),
  #              FCR = food_prov/dw) %>% 
  #       mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs), feed = as.factor(feed_types)) %>% 
  #       group_by(farm_ID, feed) %>% 
  #       mutate(prod_days = days - min(days)+1) %>% 
  #       ungroup()
  #   },
  #   pattern = cross(map(farm_IDs, farm_temp, farm_times), map(feed_types, feed_params_protein, feed_params_carbs, feed_params_lipids))
  # ),
  
  tar_target(
    farm_harvest_size,
    command = {
      df <- fish_growth(
        pop_params = pop_params,
        species_params = species_params,
        water_temp = farm_temp,
        feed_params = list(
          Proteins = feed_params_protein[[1]],
          Carbohydrates = feed_params_carbs[[1]],
          Lipids = feed_params_lipids[[1]]
        ),
        times = farm_times,
        init_weight = pop_params["meanW"],
        ingmax = pop_params["meanImax"]
      )
      df <- df %>% 
        as.data.frame() %>% remove_rownames() %>% 
        mutate(days = as.integer(days), farm_ID = as.integer(farm_IDs)) %>% 
        group_by(farm_ID) %>% 
        mutate(prod_days = days - min(days)+1) %>% 
        ungroup()
      df[df$prod_days == max(df$prod_days), ]
    },
    pattern = cross(map(farm_IDs, farm_temp, farm_times))
  ),
  
  ## Senstivities ---------------------------------------------------------------------------------------------------
  tar_target(sens_all_params, 
             command = {
               vc <- c(species_params, pop_params)
               nms <- names(vc)[!names(vc) %in% c('betaprot', 'betalip', 'betacarb', 'fcr', 'CS', 'nruns', 'deltaW', 'deltaImax', 'mortmyt', 'Wlb')]
               vc <- vc[!names(vc) %in% c('betaprot', 'betalip', 'betacarb', 'fcr', 'CS', 'nruns', 'deltaW', 'deltaImax', 'mortmyt', 'Wlb')]
               names(vc) <- nms
               vc
             }),
  tar_target(sens_params_names, names(sens_all_params)),
  tar_target(factors, c(0.9, 1, 1.1)),
  
  tar_target(
    name = sens_adjusted_params,
    command = {
      vc <- sens_all_params
      vc[sens_params_names] <- vc[sens_params_names] * factors
      vc
    },
    pattern = cross(sens_params_names, factors),
    iteration = "list"
  ),
  tar_target(
    sens_individual,
    command = {
      sens <- fish_growth(
        pop_params = sens_adjusted_params,
        species_params = sens_adjusted_params,
        water_temp = farm_temp,
        feed_params = list(
          Proteins = feed_params_protein[[1]],
          Carbohydrates = feed_params_carbs[[1]],
          Lipids = feed_params_lipids[[1]]
        ),
        times = farm_times,
        init_weight = pop_params["meanW"],
        ingmax = pop_params["meanImax"]
      ) %>% as.data.frame() %>% remove_rownames() %>% 
        mutate(adj_param = sens_params_names,
               factor = factors)
      data.frame(
        weight = last(sens$weight[!is.na(sens$weight)]),
        dw = mean(sens$dw, na.rm = T),
        P_excr = sum(sens$P_excr, na.rm = T),
        L_excr = sum(sens$L_excr, na.rm = T),
        C_excr = sum(sens$C_excr, na.rm = T),
        P_uneat = sum(sens$P_uneat, na.rm = T),
        L_uneat = sum(sens$L_uneat, na.rm = T),
        C_uneat = sum(sens$C_uneat, na.rm = T),
        food_prov = sum(sens$food_prov, na.rm = T),
        rel_feeding = mean(sens$rel_feeding, na.rm = T),
        ing_act = sum(sens$ing_act, na.rm = T),
        anab = sum(sens$anab, na.rm = T),
        catab = sum(sens$catab, na.rm = T),
        O2 = sum(sens$O2, na.rm = T),
        NH4 = sum(sens$NH4, na.rm = T),
        adj_param = unique(sens$adj_param),
        factor = unique(sens$factor)
      )
    },
    pattern = cross(map(farm_IDs, farm_temp, farm_times), map(sens_adjusted_params, cross(sens_params_names, factors)))
  )
)


