# Load required packages
suppressMessages(suppressWarnings({
  library(future)
  library(furrr)
  library(arrow)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(terra)
  library(magrittr)
  library(stringr)
  library(ggplot2)
  library(matrixStats)
  library(tibble)
  library(readxl)
  library(tictoc)
}))

conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = TRUE)

# Source model functions
source("00_model_functions.R")

# Set up parallel processing
plan(multisession, workers = future::availableCores() - 1)

# Define paths
gendata_path <- file.path("data", "_general_data")
this_species <- "atlantic_salmon"
species_path <- file.path("data", this_species)

# Load farm harvest size data (generated from mean individual modelling)
farm_harvest_size <- tar_read(farm_harvest_size, store = "04_targets_individual") %>%
  select(c(farm_ID, weight))

# Load and prepare farm data
farm_data <- read_parquet(file.path(gendata_path, "SST", "farm_SST_extracted.parquet")) %>% 
  select(data_year, tonnes_per_farm, daily_mort_rate, farm_id, day, temp_c) %>% 
  mutate(day = str_split_i(day, "day_", 2) %>% as.numeric())
farms_to_omit <- qs::qread(sprintf(file.path(gendata_path, "farm_locations", "%s_farms_to_omit.qs"), this_species))

# Get farm IDs
farm_IDs <- farm_data %>%
  filter(!farm_id %in% farms_to_omit) %>%
  distinct(farm_id) %>%
  slice_sample(n = 10) %>%
  pull(farm_id)

# Define time parameters
times_N <- c("t_start" = 121, "t_end" = 121+547, "dt" = 1)
times_S <- c("t_start" = 274, "t_end" = 274+547, "dt" = 1)

# Load and prepare farm coordinates
farm_coords <- read_parquet(file.path(gendata_path, "farm_locations", "farm_coords.parquet")) %>%
  mutate(t_start = case_when(lat > 0 ~ times_N['t_start'], TRUE ~ times_S['t_start']),
         t_end = case_when(lat > 0 ~ times_N['t_end'], TRUE ~ times_S['t_end']))

# Process farm time series data in parallel
tic()
farm_ts_data <- future_map(farm_IDs, function(id) {
  farm_data %>%
    filter(farm_id == id) %>%
    select(c(farm_id, day, temp_c))
}, .options = furrr_options(seed = TRUE))
toc() # 5.17 sec elapsed (10 farm_IDs)

# Load species parameters
species_params <- readxl::read_excel(
  file.path(species_path, "params", "Species.xlsx"),
  sheet = "Atlantic salmon"
) %>%
  {setNames(.$Value[!is.na(.$Value)], .$Quantity[!is.na(.$Value)])}

# Load population parameters
pop_params <- readxl::read_excel(file.path(species_path, "params", "Population.xlsx")) %>%
  {setNames(.$Value[!is.na(.$Value)], .$Quantity[!is.na(.$Value)])}

# Load feed parameters
feed_types <- c("reference", "past", "future")
tic()
feed_profiles <- future_map(feed_types, function(ft) {
  df <- readxl::read_excel(
    file.path(gendata_path, "diets", "feed_profiles_and_composition_Atlantic_salmon.xlsx"),
    sheet = ft
  )
  
  list(
    proteins = df %>%
      select(ingredient, proportion, contains("protein")) %>%
      select(-contains("feed")) %>%
      rename(macro = ing_protein, digest = ing_protein_digestibility),
    
    carbs = df %>%
      select(ingredient, proportion, contains("carb")) %>%
      select(-contains("feed")) %>%
      rename(macro = ing_carb, digest = ing_carb_digestibility),
    
    lipids = df %>%
      select(ingredient, proportion, contains("lipid")) %>%
      select(-contains("feed")) %>%
      rename(macro = ing_lipid, digest = ing_lipid_digestibility)
  )
}, .options = furrr_options(seed = TRUE))
toc() # 0.5 sec elapsed

# Run main farm growth calculations in parallel
tic()
main_farm_growth <- future_map(farm_IDs, function(id) {
  farm_temp <- farm_ts_data[[which(farm_IDs == id)]]$temp_c
  farm_times <- c(
    farm_coords$t_start[farm_coords$farm_ID == id],
    farm_coords$t_end[farm_coords$farm_ID == id],
    dt = 1
  )
  
  N_population <- generate_pop(
    harvest_n = (farm_data$tonnes_per_farm[farm_data$farm_id == id][1]*10^6)/
      farm_harvest_size$weight[farm_harvest_size$farm_ID == id],
    mort = pop_params['mortmyt'],
    times = farm_times
  )
  
  result <- future_map(feed_types, function(ft) {
    feed_idx <- which(feed_types == ft)
    farm_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = farm_temp,
      feed_params = list(
        Proteins = feed_profiles[[feed_idx]]$proteins,
        Carbohydrates = feed_profiles[[feed_idx]]$carbs,
        Lipids = feed_profiles[[feed_idx]]$lipids
      ),
      times = farm_times,
      N_pop = N_population/2,
      nruns = pop_params['nruns']
    )
  }, .options = furrr_options(seed = TRUE))
  
  message(paste0("Completed processing for farm_ID: ", id, " at ", Sys.time()))
  return(result)
}, .options = furrr_options(seed = TRUE))
toc() 

# Process results similar to the targets pipeline...
# (Additional processing code would follow here)

# Clean up parallel workers (reset to sequential processing)
plan(sequential)
