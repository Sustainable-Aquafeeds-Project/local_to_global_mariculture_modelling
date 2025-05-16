# Disable all linters for this file
# nolint start

library(arrow)
library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(magrittr)
library(furrr)
library(future)
library(tictoc)
library(ggplot2)
library(fs)
library(conflicted)
library(stringr)
library(readxl)
library(units)
library(qs)
library(here)
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

# Configuration ---------------------------------------------------------------------------------------------------
# Set up parallel processing
# plan(multisession, workers = parallelly::availableCores()-1)

# Basic configuration
overwrite <- F
this_species <- "atlantic_salmon"

# Directory structure
output_path <- here() %>% file.path("outputs")
gendata_path <- here() %>% file.path("data", "_general_data")
species_path <- here() %>% file.path("data", this_species)

# Input paths
input_farm_coords_path <- file.path(gendata_path, "farm_locations")
input_farm_sst_path <- file.path(gendata_path, "SST")
input_species_param_path <- file.path(species_path, "params")
input_feed_profile_path <- file.path(gendata_path, "diets")

# Output paths
output_farm_data_path <- file.path(output_path, "farm_data")
output_species_data_path <- file.path(output_path, "species_data")
output_sens_data_path <- file.path(output_path, "sensitivity_data")

# Create output directories
dir_create(c(output_farm_data_path, 
             output_species_data_path, 
             output_sens_data_path
))

# Filenames
farm_coords_file <- file.path(output_farm_data_path, "farm_coords.qs")
farm_ts_data_file <- file.path(output_farm_data_path, "farm_ts_data.qs")
species_params_file <- file.path(output_species_data_path, "species_params.qs")
sens_params_file <- file.path(output_species_data_path, "sens_params.qs")
pop_params_file <- file.path(output_species_data_path, "pop_params.qs")
feed_params_file <- file.path(output_species_data_path, "feed_params.qs")
farm_harvest_file <- file.path(output_farm_data_path, "farm_harvest_size.qs")

# Print configuration summary
message(sprintf("Configuration:\n- Species: %s\n- Parallel workers: %d\n- Overwrite existing: %s\n- Output path: %s",
               this_species,
               parallelly::availableCores()-2,
               ifelse(overwrite, "yes", "no"),
               output_path))

# Helper functions ------------------------------------------------------------------------------------------------
source("00_model_functions.R")
file_exists_skip <- function(filepath, section_name, fn) {
  if (file.exists(filepath) && !overwrite) {
    message(sprintf("Skipping %s - file exists", section_name))
    return(qs::qread(filepath))
  } else {
    tic()
    result <- fn()
    qs::qsave(result, filepath)
    message(sprintf("Saved %s - %s", section_name, toc(quiet = T)$callback_msg))
    return(result)
  }
}

# Farm coordinates ------------------------------------------------------------------------------------------------
farm_coords <- file_exists_skip(farm_coords_file, "farm coordinates", function() {
    times_N <- c("t_start" = 121, "t_end" = 121+547, "dt" = 1)
    times_S <- c("t_start" = 274, "t_end" = 274+547, "dt" = 1)
    
    read_parquet(file.path(input_farm_coords_path, "farm_coords.parquet")) %>% 
      mutate(t_start = case_when(lat > 0 ~ times_N['t_start'], TRUE ~ times_S['t_start']), 
             t_end = case_when(lat > 0 ~ times_N['t_end'], TRUE ~ times_S['t_end']),
             t_start = unname(t_start),
             t_end = unname(t_end)) 
  })

# Farm data -------------------------------------------------------------------------------------------------------
farm_ts_data <- file_exists_skip(farm_ts_data_file, "farm ts data", function() {
    farms_to_omit <- qs::qread(sprintf(file.path(input_farm_coords_path, "%s_farms_to_omit.qs"), this_species))
    
    farm_SST_data <- read_parquet(file.path(input_farm_sst_path, "farm_SST_extracted.parquet"))
    farm_IDs <- farm_SST_data %>%
      filter(!farm_id %in% farms_to_omit) %>%
      distinct(farm_id) %>%
      pull(farm_id)
    
    farm_SST_data %>%
      rename(farm_ID = farm_id) %>% 
      select(c(farm_ID, day, temp_c)) %>%
      mutate(day = str_split_i(day, "day_", 2) %>% as.numeric())
  })

# Species parameters ----------------------------------------------------------------------------------------------
species_params <- file_exists_skip(species_params_file, "species parameters", function() {
    df <- readxl::read_excel(file.path(input_species_param_path, "Species.xlsx"), sheet = "Atlantic salmon")
    vc <- df$Value
    names(vc) <- df$Quantity
    vc[!is.na(vc)]
  })

# Population parameters -------------------------------------------------------------------------------------------
pop_params <- file_exists_skip(pop_params_file, "population parameters", function() {
    df <- readxl::read_excel(file.path(input_species_param_path, "Population.xlsx"))
    vc <- df$Value
    names(vc) <- df$Quantity
    vc[!is.na(vc)]
  })

# Feed parameters -------------------------------------------------------------------------------------------------
feed_params <- file_exists_skip(feed_params_file, "feed parameters", function() {
    feed_types <- c("reference", "past", "future")
    feed_profile_file <- file.path(input_feed_profile_path, "feed_profiles_and_composition_Atlantic_salmon.xlsx")
    
    future_map(feed_types, function(ft) {
      df <- readxl::read_excel(feed_profile_file, sheet = ft)
      list(
        Proteins = df %>% select(ingredient, proportion, contains("protein")) %>%
          select(-contains("feed")) %>%
          rename(macro = ing_protein, digest = ing_protein_digestibility),
        Carbohydrates = df %>% select(ingredient, proportion, contains("carb")) %>%
          select(-contains("feed")) %>%
          rename(macro = ing_carb, digest = ing_carb_digestibility),
        Lipids = df %>% select(ingredient, proportion, contains("lipid")) %>%
          select(-contains("feed")) %>%
          rename(macro = ing_lipid, digest = ing_lipid_digestibility)
      )
    }) %>% setNames(feed_types)
  })

# Farm harvest size -----------------------------------------------------------------------------------------------
farm_harvest_size <- file_exists_skip(farm_harvest_file, "farm harvest size", function() {
    farm_IDs <- farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID)
    
    future_map_dfr(farm_IDs, function(farm_id) {
      farm_times <- c(
        t_start = farm_coords$t_start[farm_coords$farm_ID == farm_id],
        t_end = farm_coords$t_end[farm_coords$farm_ID == farm_id],
        dt = 1
      )
      
      farm_temp <- farm_ts_data %>%
        filter(farm_ID == farm_id) %>%
        filter(day >= farm_times['t_start'], day <= farm_times['t_end']) %>%
        pull(temp_c)
      
      df <- fish_growth(
        pop_params = pop_params,
        species_params = species_params,
        water_temp = farm_temp,
        feed_params = feed_params[["reference"]],
        times = farm_times,
        init_weight = pop_params["meanW"],
        ingmax = pop_params["meanImax"]
      )
      
      df %>%
        as.data.frame() %>%
        slice_tail(n = 1) %>%
        select(weight) %>%
        mutate(farm_ID = farm_id,
              weight = units::set_units(weight, "g"))
    }, .progress = TRUE)
  })

# Sensitivities ---------------------------------------------------------------------------------------------------
sens_all_params <- file_exists_skip(sens_params_file, "sensitivity parameters", function() {
    vc <- c(species_params, pop_params)
    omit_names <- c('betaprot', 'betalip', 'betacarb', 'fcr', 'CS', 'nruns')
    nms <- names(vc)[!names(vc) %in% omit_names]
    vc <- vc[!names(vc) %in% omit_names]
    names(vc) <- nms
    vc
  })

reference_feed <- feed_params[["reference"]]
factors <- c(0.9, 1, 1.1)
sens_params_names <- names(sens_all_params)
farm_IDs <- farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID) %>% sample(10)

rm(list = grep("tar_", ls(), value = TRUE), envir = .GlobalEnv)
targets::tar_make(
  store = "_targets_sensitivities",
  script = "_targets_sensitivities.R",
  reporter = "balanced",
  callr_function = NULL
)

overwrite <- T
sens_results <- tar_read(tar_sens_results, store = "_targets_sensitivities") %>% 
  mutate(measure = as.factor(measure),
         adj_param = factor(adj_param, levels = sens_params_names))
sens_measures <- levels(sens_results$measure)
sens_results_files <- file.path(output_sens_data_path, paste0("sens_results_", sens_measures, ".qs"))
sens_results_figfiles <- file.path(output_sens_data_path, paste0("sens_plot_", sens_measures, ".qs"))
for (sm in seq_along(sens_measures)) {
  file_exists_skip(sens_results_files[sm], paste("sensitivity files", sm, "of", length(sens_results_files)), function() {
    sens_results %>% 
      filter(measure == sens_measures[sm])
  })
  file_exists_skip(sens_results_figfiles[sm], paste("sensitivity plots", sm, "of", length(sens_results_figfiles)), function() {
    sens_results %>% 
      filter(measure == sens_measures[sm]) %>% 
      ggplot(aes(x = adj_param, y = mean_sens, ymin = mean_sens-sd_sens, ymax = mean_sens+sd_sens)) +
      geom_col(fill = "salmon", alpha = 0.65, colour = "black") +
      geom_errorbar(width = 0.5) +
      coord_flip() +
      theme_classic()
  })
}

# Finish up -------------------------------------------------------------------------------------------------------
plan(sequential)
message("All done!")

# nolint end
