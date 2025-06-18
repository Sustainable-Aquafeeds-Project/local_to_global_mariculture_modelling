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

here("00_model_functions.R") %>% source()
here("00_dirs.R") %>% source()

# Configuration ---------------------------------------------------------------------------------------------------
# Set up parallel processing
# plan(multisession, workers = parallelly::availableCores()-2)
# plan(multisession, workers = works)

# Filenames
farm_geometry_file <- file.path(output_farm_data_path, "farm_geometry.qs")
sens_params_file <- file.path(output_species_data_path, "sens_params.qs")

# Print configuration summary
message(sprintf("Configuration:\n- Species: %s\n- Overwrite existing: %s\n- Output path: %s",
               this_species,
               ifelse(overwrite, "yes", "no"),
               output_path))

# Main farm growth ------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_farmruns")
rm(list = grep("tar_", ls(), value = TRUE), envir = .GlobalEnv)
targets::tar_make(seconds_meta_append = 90)

tar_farm_IDs <- targets::tar_read(tar_farm_IDs)
stat_names <- targets::tar_read(stat_names)
farmrun_reference_files <- file.path(output_growth_data_path, sprintf("farmrun_reference_%s.qs", fixnum(tar_farm_IDs,4)))
farmrun_past_files <- file.path(output_growth_data_path, sprintf("farmrun_past_%s.qs", fixnum(tar_farm_IDs,4)))
farmrun_future_files <- file.path(output_growth_data_path, sprintf("farmrun_future_%s.qs", fixnum(tar_farm_IDs,4)))
farmrun_comparisons_files <- file.path(output_cohorts_data_path, sprintf("farmrun_comparisons_%s.qs", fixnum(tar_farm_IDs,4)))

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
  
  # Save the final cohort feed comparison data
  targets::tar_read(tar_farmrun_comparisons, branches = f) %>% 
    setNames(stat_names) %>% 
    qs::qsave(farmrun_comparisons_files[f])
}

# Use the intermediate growth data to save aggregate files
aggregate_comparison_files <- file.path(output_cohorts_data_path, str_c("allfarms_comparisons_", stat_names, ".qs"))

for (s in seq_along(stat_names)) {
  purrr::map(farmrun_comparisons_files, function(fnm) {
    qread(fnm)[[stat_names[s]]]
  }) %>% 
    bind_rows() %>% 
    qsave(aggregate_comparison_files[s])
}

# nolint end
