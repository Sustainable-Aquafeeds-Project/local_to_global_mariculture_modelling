# Setup -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
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
})

# functions
source("00_model_functions.R")

# species paths
this_species <- "atlantic_salmon"
this_path <- file.path("data", this_species)
fig_path <- file.path(this_path, "figures")
prod_path <- file.path(this_path, "data_products")

# STEP 1: Temperature forcings ------------------------------------------------------------------------------------
# You don't need to do this again unless you make changes - all the important stuff is saved
# Prep temperature forcings for each farm site
# production_cycle <- read.csv("data/_general_data/production_cycles/production_cycle.csv") %>% 
#   filter(species == this_species) %>% 
#   rename(production_cycle_length = days) %>% 
#   pull(production_cycle_length)
production_cycle <- 1100

farms <-  qread("data/_general_data/farm_locations/locations_w_species_fao_area_stocking.qs") %>% 
  filter(model_name == this_species) %>% 
  select(-row_num) %>% 
  mutate(farm_id = row_number())

hemi <- cbind(farms$farm_id, sf::st_coordinates(farms$geometry)) %>% 
  as.data.frame() %>% rename(farm_ID = V1, lon = X, lat = Y) %>% 
  write_parquet("data/_general_data/farm_locations/farm_coords.parquet")

day_number <- seq(1:production_cycle)

temp_data <- purrr::map_dfc(.x = day_number, .f = function(day_number){
  rast_day_number <- if_else(day_number <= 365, true = day_number, false = day_number-365)
  rast_day_number <- if_else(rast_day_number <= 365, true = rast_day_number, false = rast_day_number-365)
  rast_day_number <- if_else(rast_day_number <= 365, true = rast_day_number, false = rast_day_number-365)
  message("Getting temperature data for all sites for ", this_species,  " - day ", day_number)
  
  sst_test <- terra::rast(sprintf("data/_general_data/SST/SST_gf_rasters/sst_nasa_mur_L4_0.25_mean2010-2019_day_%s.tif", rast_day_number))
  
  terra::extract(sst_test, farms) %>%
    mutate(day = paste0("day_", day_number)) %>%
    pivot_wider(names_from = "day", values_from = "focal_mean") %>%
    select(-ID)
}) %>%
  mutate(farm_id = row_number())
# If you want the sf object it's here!

farms_w_temp_df <- farms %>%
  left_join(temp_data, by = c("farm_id" = "farm_id")) %>%
  pivot_longer(names_to = "day", values_to = "temp_c", cols = starts_with("day_"))

# Check which farms have missing temp data
(
  missing_temp_farms <- farms_w_temp_df %>% 
    filter(temp_c %>% is.na()) %>% 
    group_by(farm_id) %>% 
    reframe(num_missing = n())
)

# How far apart in the sequence are the farms? If the previous is complete we should be able to use the one before in the same country
diff(missing_temp_farms$farm_id)

# Make the farm list
farm_list <- farms_w_temp_df %>%
  group_by(farm_id) %>% 
  group_split()

# Loop through and assigned temp of farms missing temp data, to the farm adjacent (the nearest complete index before)
for(i in 1:length(farm_list)){
  message("Checking temp data for ", unique(farm_list[[i]]$farm_id)) 
  if(unique(is.na(farm_list[[i]]$temp_c))){ #if temp data is NA see below
    cat("Is the previous farm index the same country?")
    if(unique(farm_list[[i-1]]$country) == unique(farm_list[[i]]$country)){
      if(!unique(is.na(farm_list[[i-1]]$temp_c))){ # if the farm index before is NOT NA, use that.
        farm_list[[i]]$temp_c <- farm_list[[i-1]]$temp_c
      } else {
        farm_list[[i]]$temp_c <- farm_list[[i-2]]$temp_c.  #else use the farm index 2 before (the missing_farm_
      }
    } else {stop("Previous country index not the same")} #if the previous country is not the same country stop the loop
  }
}

# Check again - looks good - no values.
bind_rows(farm_list) %>%  filter(temp_c %>% is.na()) %>% pull(farm_id) %>% unique()

# Save the new locations data 
farms_w_temp_df <- bind_rows(farm_list)

# With geometry, for plotting
qsave(x = farms_w_temp_df, 
      file = sprintf("data/_general_data/farm_locations/%s_locations_w_temps.qs", this_species))

# Without geometry, for targets
sf::st_drop_geometry(farms_w_temp_df) %>%
  write_parquet("data/_general_data/SST/farm_SST_extracted.parquet")

# Get the mean temps for each farm - this is needed to check which annual temps the fish can deal with
mean_farm_temp <- farm_list %>% 
  map_df(.f = function(x){
    data.frame(farm_id = unique(x$farm_id), 
               mean_temp = mean(x$temp_c),
               country = unique(x$country),
               volume = unique(x$tonnes_per_farm))
  })

farms_to_omit <- mean_farm_temp %>% 
  filter(mean_temp <= 6) %>% 
  pull(farm_id)

qsave(x = farms_to_omit, 
      file = sprintf("data/_general_data/farm_locations/%s_farms_to_omit.qs", this_species))


# STEP 2 - Run model ----------------------------------------------------------------------------------------------
## Example individuals --------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_individual")
# tar_visnetwork()
tar_make(reporter = "balanced", seconds_meta_append = 300)
# tar_prune()

farm_IDs <- tar_read(farm_IDs)
sens_params_names <- tar_read(sens_params_names)
patt <- data.frame(
  farm_ID = rep(farm_IDs, each = length(sens_params_names)*3),
  param = rep(rep(1:length(sens_params_names), each = 3), times = length(farm_IDs)),
  factor = rep(c(0.9,1,1.1), times = length(sens_params_names)*length(farm_IDs))
)
patt$br <- 1:nrow(patt)

wt_ls <- dw_ls <- excr_ls <- uneat_ls <- list()
for (p in 1:length(sens_params_names)) {
  br <- patt[patt$param == p, ]
  sens <- tar_read(sens_individual, branches = br$br)
  sens$farm_ID <- br$farm_ID
  
  wt_ls[[p]] <- sens %>% 
    select(weight, farm_ID, adj_param, factor) %>% 
    pivot_wider(names_from = factor, names_prefix = "p", values_from = weight, id_cols = c(adj_param, farm_ID)) %>% 
    mutate(sens = (p1.1 - p0.9)/(0.2*p1)) %>% 
    group_by(adj_param) %>% 
    reframe(sd = sd(sens), sens = mean(sens))
  
  dw_ls[[p]] <- sens %>% 
    select(dw, farm_ID, adj_param, factor) %>% 
    pivot_wider(names_from = factor, names_prefix = "p", values_from = dw) %>% 
    mutate(sens = (p1.1 - p0.9)/(0.2*p1)) %>% 
    group_by(adj_param) %>% 
    reframe(sd = sd(sens), sens = mean(sens))
  
  excr_ls[[p]] <- sens %>% 
    mutate(excr = P_excr + L_excr + C_excr) %>% 
    select(excr, farm_ID, adj_param, factor) %>% 
    pivot_wider(names_from = factor, names_prefix = "p", values_from = excr) %>% 
    mutate(sens = (p1.1 - p0.9)/(0.2*p1)) %>% 
    group_by(adj_param) %>% 
    reframe(sd = sd(sens), sens = mean(sens))
  
  uneat_ls[[p]] <- sens %>% 
    mutate(uneat = P_uneat + L_uneat + C_uneat) %>% 
    select(uneat, farm_ID, adj_param, factor) %>% 
    pivot_wider(names_from = factor, names_prefix = "p", values_from = uneat) %>% 
    mutate(sens = (p1.1 - p0.9)/(0.2*p1)) %>% 
    group_by(adj_param) %>% 
    reframe(sd = sd(sens), sens = mean(sens))
}
wt_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "sensitivities", "weight_parameter_sensitivity.parquet"))
dw_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "sensitivities", "dw_parameter_sensitivity.parquet"))
excr_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "sensitivities", "excreted_parameter_sensitivity.parquet"))
uneat_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "sensitivities", "uneaten_parameter_sensitivity.parquet"))

## Farm growth ----------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_farm")
# tar_visnetwork()
tar_validate()
# tar_prune()
upd <- tar_outdated(branches = F)
tar_make(reporter = "balanced", seconds_meta_append = 300)

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
