# Setup -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(terra)
  library(qs)
  library(here)
  library(arrow)
  library(units)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(conflicted)
  conflicts_prefer(dplyr::select(), dplyr::filter(), .quiet = T)
})

# This script aims to take the outputs from the model and answer specific questions 
# Paths & globals -------------------------------------------------------------------------------------------------
source("00_model_functions.R")

coho_path <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohorts")

# * Main point is the difference between feeds. 
# * This should be expressed in t/t salmon, so farm biomass is also important.
# * Also interesting to see if difference between feeds varies geographically (correlates with mean/median/max temperature?)

# Biomass produced ------------------------------------------------------------------------------------------------
# How accurately is the model producing the correct biomass production levels?
# Here, "cohorts" is going to become "farms", or more accurately, farm outputs per 2-year cycle with overlaps (cohort 2) 

farm_info_ts <- file.path("data", "_general_data", "farm_locations", "atlantic_salmon_locations_w_temps.qs") %>% 
  qread() %>% 
  filter(day == "day_1") %>% 
  mutate(country = as.factor(country),
         species_group = as.factor(species_group),
         data_type = as.factor(data_type),
         data_type_2 = as.factor(data_type_2),
         model_name = as.factor(model_name),
         F_CODE = as.factor(F_CODE),
         day = str_split_i(day, "_", 2) %>% as.integer())
hist(farm_info$tonnes_per_farm)


coho_biom <- file.path(coho_path) %>% 
  list.files(full.names = T) %>% 
  str_subset("biomass")



# Difference between feeds ----------------------------------------------------------------------------------------




