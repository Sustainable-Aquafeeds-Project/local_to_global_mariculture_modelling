library(targets)
library(dplyr)
library(arrow)
library(tictoc)
library(furrr)
library(future)

source("00_model_functions.R")

farm_IDs <- tar_read(farm_IDs)
feed_types <- tar_read(feed_types)

data_path_farm <- file.path(prod_path, "model_outputs_farm")
data_path_cohort <- file.path(prod_path, "model_outputs_cohort")
dest_path_farm <- file.path(prod_path, "processed_farm")
if (!dir.exists(dest_path_farm)) {dir.create(dest_path_farm)}
dest_path_cohort <- file.path(prod_path, "processed_cohort")
if (!dir.exists(dest_path_cohort)) {dir.create(dest_path_cohort)}

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 10)

# Farm processing - not currently in use
# farm_stats <- c("weight_stat", "biomass_stat", "dw_stat", "SGR_stat", "E_somat_stat", "P_excr_stat", "L_excr_stat", 
#                 "C_excr_stat", "P_uneat_stat", "L_uneat_stat", "C_uneat_stat", "ing_act_stat", "anab_stat", 
#                 "catab_stat", "NH4_stat", "O2_stat", "food_prov_stat", "rel_feeding_stat", "T_response_stat")
# # List all the files output from targets processing
# farm_files <- list()
# for (s in seq_along(farm_stats)) {
#   farm_files[[s]] <- list.files(data_path_farm, full.names = T) %>% 
#     str_subset(farm_stats[s])
# }
# # Processing here?
# rm(farm_files)

# Cohorts processing ---------------------------------------------------------------------------------------------------
cohort_files <- list.files(data_path_cohort, full.names = T)

## Uneaten feed --------------------------------------------------------------------------------------------------------
# Function to combine uneaten feed
process_uneat <- function(fid) {
  par_name <- file.path(dest_path_cohort, paste0("total_uneat_stat_", fixnum(farm_IDs[fid]), ".parquet"))
  if (overwrite == T | !file.exists(par_name)) {
    fnms <- cohort_files %>% 
      str_subset(fixnum(farm_IDs[fid])) %>% 
      str_subset("uneat")
    
    df <- list(
      fnms %>% str_subset("C_uneat") %>% read_parquet() %>% mutate(stat = "C_uneat"),
      fnms %>% str_subset("L_uneat") %>% read_parquet() %>% mutate(stat = "L_uneat"),
      fnms %>% str_subset("P_uneat") %>% read_parquet() %>% mutate(stat = "P_uneat")
    ) %>% 
      bind_rows() %>% 
      dplyr::filter(yday != max(yday)) %>% 
      mutate(stat = as.factor(stat),
             cohort = as.factor(cohort)) %>% 
      write_parquet(par_name)
  }
}

# Process uneaten stats in parallel
future_map(seq_along(farm_IDs), process_uneat, .progress = TRUE)

## Excreted faeces -----------------------------------------------------------------------------------------------------
# Function to combine excreted
process_excr <- function(fid) {
  par_name <- file.path(dest_path_cohort, paste0("total_excr_stat_", fixnum(farm_IDs[fid]), ".parquet"))
  if (overwrite == T | !file.exists(par_name)) {
    fnms <- cohort_files %>% 
      str_subset(fixnum(farm_IDs[fid])) %>% 
      str_subset("excr")
    
    df <- list(
      fnms %>% str_subset("C_excr") %>% read_parquet() %>% mutate(stat = "C_excr"),
      fnms %>% str_subset("L_excr") %>% read_parquet() %>% mutate(stat = "L_excr"),
      fnms %>% str_subset("P_excr") %>% read_parquet() %>% mutate(stat = "P_excr")
    ) %>% 
      bind_rows() %>% 
      dplyr::filter(yday != max(yday)) %>% 
      mutate(stat = as.factor(stat),
             cohort = as.factor(cohort)) %>% 
      write_parquet(par_name)
  }
}

# Process excreted stats in parallel
future_map(seq_along(farm_IDs), process_excr, .progress = TRUE)

# Clean up parallel processing
plan(sequential)
