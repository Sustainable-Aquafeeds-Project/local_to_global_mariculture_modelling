# Setup -----------------------------------------------------------------------------------------------------------
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(terra)
library(qs)
library(here)
library(sf)
library(purrr)
library(furrr)
library(targets)
library(future)
library(arrow)
library(readxl)
library(units)
library(tictoc)
library(conflicted)
conflicts_prefer(dplyr::select(), dplyr::filter(), .quiet = T)

# functions
source("00_model_functions.R")

# species paths
this_species <- "atlantic_salmon"
this_path <- file.path("data", this_species)
fig_path <- file.path(this_path, "figures")
prod_path <- file.path(this_path, "data_products")

# STEP 1: Temperature forcings ------------------------------------------------------------------------------------
# You don't need to do this again unless you make changes - all the important stuff is already saved
source("04_extracting_temperatures.R")

# STEP 2 - Run model ----------------------------------------------------------------------------------------------
## Individuals ----------------------------------------------------------------------------------------------------
source("06.1_run_individual.R")

### Harvest size ---------------------------------------------------------------------------------------------------
tar_read(farm_harvest_size, store = "04_targets_individual") %>% 
  write_parquet(file.path("data", "atlantic_salmon", "data_products", "harvest_size.parquet"))

## Farm growth ----------------------------------------------------------------------------------------------------
source("06.2_run_farms.R")

# source("10_process_targets_outputs.R")
# Analysis --------------------------------------------------------------------------------------------------------

## Combining previously processed data ----------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_farm")

overwrite <- T
source("11_process_targets_outputs_2.R")

# Final data outputs ----------------------------------------------------------------------------------------------
excr_fnms <- list.files("data/atlantic_salmon/data_products/processed_cohort", full.names = T) %>% str_subset("total_excr")

excr_ls <- uneat_ls <- list()
for (f in seq_along(excr_fnms)) {
  excr_ls[[f]] <- excr_fnms[f] %>% arrow::read_parquet()
  uneat_ls[[f]] <- uneat_fnms[f] %>% arrow::read_parquet()
}
excr_ls %>% bind_rows() %>% 
  arrow::write_parquet("data/atlantic_salmon/data_products/total_excr_allfarms.parquet")
uneat_ls %>% bind_rows() %>% 
  arrow::write_parquet("data/atlantic_salmon/data_products/total_uneat_allfarms.parquet")

source("12_data_analysis.R")

# Benchmarking ----------------------------------------------------------------------------------------------------
# rbenchmark::benchmark(
#   df = {apportion_feed(provided = 10, ingested = 9,
#                            prop = feed_params[['Carbohydrates']]$proportion,
#                            macro = feed_params[['Carbohydrates']]$ing_carb,
#                            digestibility = feed_params[['Carbohydrates']]$ing_carb_digestibility)},
#   matrix = {apportion_feed_short(provided = 10, ingested = 9,
#                              prop = feed_params[['Carbohydrates']]$proportion,
#                              macro = feed_params[['Carbohydrates']]$ing_carb,
#                              digestibility = feed_params[['Carbohydrates']]$ing_carb_digestibility)},
#   replications = 1000
# )
