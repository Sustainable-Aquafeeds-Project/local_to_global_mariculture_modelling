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
library(fs)
library(conflicted)
library(stringr)
library(readxl)
library(units)
library(qs)
library(here)
library(progressr)
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

# Configuration ---------------------------------------------------------------------------------------------------
# Set up parallel processing
# plan(multisession, workers = parallelly::availableCores()-2)
# plan(multisession, workers = works)

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
output_growth_data_path <- file.path(output_path, "farm_growth_data")
output_cohorts_data_path <- file.path(output_path, "cohort_growth_data")
output_model_farm_path <- file.path(output_path, "all_outputs_farm")
output_model_cohort_path <- file.path(output_path, "all_outputs_cohort")
total_uneaten_path <- file.path(output_path, "total_uneaten_cohort")
total_excreted_path <- file.path(output_path, "total_excreted_cohort")

# Filenames
farm_coords_file <- file.path(output_farm_data_path, "farm_coords.qs")
farm_geometry_file <- file.path(output_farm_data_path, "farm_geometry.qs")
farm_ts_data_file <- file.path(output_farm_data_path, "farm_ts_data.qs")
species_params_file <- file.path(output_species_data_path, "species_params.qs")
sens_params_file <- file.path(output_species_data_path, "sens_params.qs")
pop_params_file <- file.path(output_species_data_path, "pop_params.qs")
feed_params_file <- file.path(output_species_data_path, "feed_params.qs")
farm_harvest_file <- file.path(output_farm_data_path, "farm_harvest_size.qs")

# Create output directories
dir_create(c(output_farm_data_path, 
             output_species_data_path, 
             output_growth_data_path,
             output_cohorts_data_path,
             output_model_farm_path,
             output_model_cohort_path,
             total_uneaten_path,
             total_excreted_path))

# Print configuration summary
message(sprintf("Configuration:\n- Species: %s\n- Overwrite existing: %s\n- Output path: %s",
               this_species,
               ifelse(overwrite, "yes", "no"),
               output_path))

# Helper functions ------------------------------------------------------------------------------------------------
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

source("00_model_functions.R")

# Load farm harvest size data -------------------------------------------------------------------------------------
# This is from the individual running
farm_harvest_size <- qs::qread(file.path(output_farm_data_path, "farm_harvest_size.qs")) %>%
  select(c(farm_ID, weight)) %>%
  mutate(weight = units::set_units(weight, "g"))

# Other previously saved data -------------------------------------------------------------------------------------
farms_to_omit <- qs::qread(sprintf(file.path(input_farm_coords_path, "%s_farms_to_omit.qs"), this_species))
farm_coords <- qs::qread(farm_coords_file) %>% filter(!farm_ID %in% farms_to_omit)
farm_ts_data <- qs::qread(farm_ts_data_file) %>% filter(!farm_ID %in% farms_to_omit)
species_params <- qs::qread(species_params_file)
pop_params <- qs::qread(pop_params_file)
feed_params <- qs::qread(feed_params_file)

# Main farm growth ------------------------------------------------------------------------------------------------
stat_names <- c("weight_stat", "dw_stat", "water_temp_stat", "T_response_stat", "P_excr_stat", "L_excr_stat",
                "C_excr_stat", "P_uneat_stat", "L_uneat_stat", "C_uneat_stat", "food_prov_stat", "food_enc_stat", 
                "rel_feeding_stat", "ing_pot_stat", "ing_act_stat", "E_assim_stat", "E_somat_stat", "anab_stat",
                "catab_stat", "O2_stat", "NH4_stat", "total_excr_stat", "total_uneat_stat", "metab_stat",
                "biomass_stat")

Sys.setenv(TAR_PROJECT = "project_farmruns")
rm(list = grep("tar_", ls(), value = TRUE), envir = .GlobalEnv)
targets::tar_make(
  reporter = "balanced",
  callr_function = NULL,
  seconds_meta_append = 90
)

tar_farm_IDs <- targets::tar_read(tar_farm_IDs)
farmrun_reference_files <- file.path(output_growth_data_path, sprintf("farmrun_reference_%s.qs", fixnum(tar_farm_IDs,4)))
farmrun_past_files <- file.path(output_growth_data_path, sprintf("farmrun_past_%s.qs", fixnum(tar_farm_IDs,4)))
farmrun_future_files <- file.path(output_growth_data_path, sprintf("farmrun_future_%s.qs", fixnum(tar_farm_IDs,4)))

# Save the intermediate farm growth data
for (f in seq_along(tar_farm_IDs)) {
  targets::tar_read(tar_farmrun_reference, branches = f) %>%
    setNames(stat_names) %>% 
    lapply(function(x) {colnames(x) <- c("t", "mean", "sd"); x}) %>%
    lapply(function(x) {cbind(x, farm_ID = tar_farm_IDs[f])})  %>%
    qs::qsave(farmrun_reference_files[f])
  targets::tar_read(tar_farmrun_past, branches = f) %>% 
    setNames(stat_names) %>% 
    lapply(function(x) {colnames(x) <- c("t", "mean", "sd"); x}) %>%
    lapply(function(x) {cbind(x, farm_ID = tar_farm_IDs[f])})  %>%
    qs::qsave(farmrun_past_files[f])
  targets::tar_read(tar_farmrun_future, branches = f) %>% 
    setNames(stat_names) %>% 
    lapply(function(x) {colnames(x) <- c("t", "mean", "sd"); x}) %>%
    lapply(function(x) {cbind(x, farm_ID = tar_farm_IDs[f])})  %>%
    qs::qsave(farmrun_future_files[f])
}

farmrun_comparisons_files <- file.path(output_cohorts_data_path, sprintf("farmrun_comparisons_%s.qs", fixnum(tar_farm_IDs,4)))
# Save the final cohort feed comparison data
for (f in seq_along(tar_farm_IDs)) {
  targets::tar_read(tar_farmrun_comparisons, branches = f) %>% 
    setNames(stat_names) %>% 
    qs::qsave(farmrun_comparisons_files[f])
}

# nolint end
