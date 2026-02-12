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
  library(conflicted)
}))

conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = TRUE)

# Source model functions
source("00_model_functions.R")

# Set up parallel processing
plan(multisession, workers = future::availableCores() - 2)

# Define paths
gendata_path <- file.path("data", "_general_data")
this_species <- "atlantic_salmon"
species_path <- file.path("data", this_species)

# Load farm harvest size data (generated from mean individual modelling)
farm_harvest_size <- targets::tar_read(farm_harvest_size, store = "04_targets_individual") %>%
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
main_farm_growth <- future_map(farm_IDs, function(id) {
  farm_times <- c(
    farm_coords$t_start[farm_coords$farm_ID == id],
    farm_coords$t_end[farm_coords$farm_ID == id],
    dt = 1
  )
  farm_temp <- farm_ts_data[[which(farm_IDs == id)]]$temp_c[farm_times['t_start']:farm_times['t_end']]
  
  N_population <- generate_pop(
    harvest_n = (farm_data$tonnes_per_farm[farm_data$farm_id == id][1]*10^6)/
      farm_harvest_size$weight[farm_harvest_size$farm_ID == id],
    mort = pop_params['mortmyt'],
    times = farm_times
  )
  
  if (!file.exists(file.path("outputs", "main_farm_growth", paste0("farm_ID_", id, ".Rdata"))) | overwrite == T) {
    tic()
    result <- list()
    for (ft in seq_along(feed_types)) {
      result[[ft]] <- farm_growth(
        pop_params = pop_params,
        species_params = species_params,
        water_temp = farm_temp,
        feed_params = list(
          Proteins = feed_profiles[[ft]]$proteins,
          Carbohydrates = feed_profiles[[ft]]$carbs,
          Lipids = feed_profiles[[ft]]$lipids
        ),
        times = farm_times,
        N_pop = N_population/2,
        nruns = 50 #pop_params['nruns']
      )
    }
    save(result, file.path("outputs", "main_farm_growth", paste0("farm_ID_", id, ".Rdata")))
    print(paste0("Completed modelling of farm_ID ", id, "; ", toc(quiet = T)$callback_msg))
  } else {
    load(file.path("outputs", "main_farm_growth", paste0("farm_ID_", id, ".Rdata")))
    print(paste0("Modelling of farm_ID ", id, " loaded from previous run."))
  }
  
  # tar_target(
  #   cohorts_biomass,
  #   command = {
  #     df0 <- cbind(matrix(121:(121+547),548,1), main_farm_growth[[2]])
  #     df1 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365, cohort = 2, prod_days = 1:nrow(df0))
  #     df2 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(V1 = V1+365+365, cohort = 3, prod_days = 1:nrow(df0))
  #     df0 <- df0 %>% as.data.frame() %>% remove_rownames() %>% mutate(cohort = 1, prod_days = 1:nrow(df0))
  #     df <- rbind(df0, df1, df2) %>% rename(days = V1, mean = V2, sd = V3) %>% relocate(prod_days, .after = days) %>% 
  #       mutate(farm_ID = farm_IDs, feed = feed_types)
  #   },
  #   pattern = map(main_farm_growth, cross(feed_types, farm_IDs)), 
  #   deployment = "main"
  # ),

}, .options = furrr_options(seed = TRUE), .progress = T)

# Process results similar to the targets pipeline...
# (Additional processing code would follow here)

# Clean up parallel workers (reset to sequential processing)
plan(sequential)

# Time testing
tic()
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
    nruns = 50 #pop_params['nruns']
  )
}, .progress = T)
t1 <- toc(quiet = T)$callback_msg

tic()
t2 <- toc(quiet = T)$callback_msg

# Population scaling ----------------------------------------------------------------------------------------------
# Do I actually need 500 runs?
Sys.setenv(TAR_PROJECT = "project_farmruns")
tar_load(farm_ID_data, branches = 500)
names(farm_ID_data) <- c("t_start", "t_end", "dt", "nruns", "prod")
tar_load(N_pop, branches = 500)
tar_load(farm_ts_data, branches = 500)
tar_load(all_params)
tar_load(feed_params, branches = 1)
names(feed_params) <- c("Proteins", "Carbohydrates", "Lipids")
days <- (farm_ID_data['t_start']:farm_ID_data['t_end'])*farm_ID_data['dt']

nruns <-seq(500, 5000, 250)
results <- purrr::map(nruns, function(nr) {
  # Generate all random values upfront
  init_weights <- rnorm(nr, mean = all_params['meanW'], sd = all_params['deltaW'])
  ingmaxes <- rnorm(nr, mean = all_params['meanImax'], sd = all_params['deltaImax'])
  
  # Run parallel simulation for individuals
  mc_results <- purrr::map2(init_weights, ingmaxes, function(init_w, ing_m) {
    mat <- fish_growth(
      pop_params = all_params,
      species_params = all_params,
      water_temp = farm_ts_data$temp_c,
      feed_params = feed_params,
      times = farm_ID_data,
      init_weight = init_w,
      ingmax = ing_m
    ) %>% unname()
  })
  
  lapply(mc_results, function(mc) {
    weight <- mc[548, 2]
    P_excr <- mc[, 6] %>% sum()
    P_uneat <- mc[, 9] %>% sum()
    ing_act <- mc[, 16] %>% sum()
    cbind(weight, P_excr, P_uneat, ing_act) %>% 
      as.data.frame()
  }) %>% 
    bind_rows() %>% 
    slice_sample(n = 500)
})

# See if there are appreciable differences
weight.stats <- data.frame(
  "nruns" = nruns,
  "weight.mean" = lapply(results, function(rs) {mean(rs$weight)}) %>% unlist(),
  "weight.sd" = lapply(results, function(rs) {sd(rs$weight)}) %>% unlist()
)
for (i in seq_along(nruns)) {
  weight.stats$p.value[i] <- t.test(results[[length(nruns)]]$weight, results[[i]]$weight)[["p.value"]]
}
ggplot(weight.stats, aes(x = nruns, y = weight.mean, colour = p.value)) +
  geom_point(size = 4)

P_excr.stats <- data.frame(
  "nruns" = nruns,
  "P_excr.mean" = lapply(results, function(rs) {mean(rs$P_excr)}) %>% unlist(),
  "P_excr.sd" = lapply(results, function(rs) {sd(rs$P_excr)}) %>% unlist()
)
for (i in seq_along(nruns)) {
  P_excr.stats$p.value[i] <- t.test(results[[length(nruns)]]$P_excr, results[[i]]$P_excr)[["p.value"]]
}
ggplot(P_excr.stats, aes(x = nruns, y = P_excr.mean, colour = p.value, size = p.value)) +
  geom_point()

