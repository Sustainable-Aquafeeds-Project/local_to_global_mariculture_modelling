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

# farm_IDs <- tar_read(farm_IDs)
# sens_params_names <- tar_read(sens_params_names)
# patt <- data.frame(
#   farm_ID = rep(farm_IDs, each = length(sens_params_names)*3),
#   param = rep(rep(1:length(sens_params_names), each = 3), times = length(farm_IDs)),
#   factor = rep(c(0.9,1,1.1), times = length(sens_params_names)*length(farm_IDs))
# )
# patt$br <- 1:nrow(patt)
# 
# wt_ls <- dw_ls <- excr_ls <- uneat_ls <- list()
# for (p in 1:length(sens_params_names)) {
#   br <- patt[patt$param == p, ]
#   sens <- tar_read(sens_individual, branches = br$br)
#   sens$farm_ID <- br$farm_ID
#   
#   wt_ls[[p]] <- sens %>% 
#     select(weight, farm_ID, adj_param, factor) %>% 
#     pivot_wider(names_from = factor, names_prefix = "p", values_from = weight, id_cols = c(adj_param, farm_ID)) %>% 
#     mutate(sens = (p1.1 - p0.9)/(0.2*p1)) %>% 
#     group_by(adj_param) %>% 
#     reframe(sd = sd(sens), sens = mean(sens))
#   
#   dw_ls[[p]] <- sens %>% 
#     select(dw, farm_ID, adj_param, factor) %>% 
#     pivot_wider(names_from = factor, names_prefix = "p", values_from = dw) %>% 
#     mutate(sens = (p1.1 - p0.9)/(0.2*p1)) %>% 
#     group_by(adj_param) %>% 
#     reframe(sd = sd(sens), sens = mean(sens))
#   
#   excr_ls[[p]] <- sens %>% 
#     mutate(excr = P_excr + L_excr + C_excr) %>% 
#     select(excr, farm_ID, adj_param, factor) %>% 
#     pivot_wider(names_from = factor, names_prefix = "p", values_from = excr) %>% 
#     mutate(sens = (p1.1 - p0.9)/(0.2*p1)) %>% 
#     group_by(adj_param) %>% 
#     reframe(sd = sd(sens), sens = mean(sens))
#   
#   uneat_ls[[p]] <- sens %>% 
#     mutate(uneat = P_uneat + L_uneat + C_uneat) %>% 
#     select(uneat, farm_ID, adj_param, factor) %>% 
#     pivot_wider(names_from = factor, names_prefix = "p", values_from = uneat) %>% 
#     mutate(sens = (p1.1 - p0.9)/(0.2*p1)) %>% 
#     group_by(adj_param) %>% 
#     reframe(sd = sd(sens), sens = mean(sens))
# }
# wt_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "sensitivities", "weight_parameter_sensitivity.parquet"))
# dw_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "sensitivities", "dw_parameter_sensitivity.parquet"))
# excr_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "sensitivities", "excreted_parameter_sensitivity.parquet"))
# uneat_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "sensitivities", "uneaten_parameter_sensitivity.parquet"))

## Farm growth ----------------------------------------------------------------------------------------------------
source("06.2_run_farms.R")

overwrite <- T
source("10_process_targets_outputs.R")

# Analysis --------------------------------------------------------------------------------------------------------
## Harvest size ---------------------------------------------------------------------------------------------------
pdata <- tar_read(farm_harvest_size, store = "04_targets_individual") %>% 
  write_parquet(file.path("data", "atlantic_salmon", "data_products", "harvest_size.parquet"))
rm(pdata)

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
