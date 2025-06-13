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

# paths
source("00_dirs.R")

# STEP 1: Temperature forcings ------------------------------------------------------------------------------------
# You don't need to do this again unless you make changes - all the important stuff is already saved
source("04_extracting_temperatures.R")

# STEP 2: Formulate feeds -----------------------------------------------------------------------------------------
source("05_formulating_feeds.R")

# STEP 3: Run model -----------------------------------------------------------------------------------------------
## Individuals ----------------------------------------------------------------------------------------------------
overwrite <- F
source("06.1_run_individual.R")

## Farm growth ----------------------------------------------------------------------------------------------------
overwrite <- F
source("06.2_run_farms.R")

# STEP 3: Analysis ------------------------------------------------------------------------------------------------
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
