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
