library(targets)
library(dplyr)
library(arrow)
library(tictoc)
library(furrr)
library(future)

source("00_model_functions.R")

farm_IDs <- tar_read(farm_IDs)
feed_types <- tar_read(feed_types)

dest_path_farm <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_farm")
if (!dir.exists(dest_path_farm)) {dir.create(dest_path_farm)}
dest_path_cohort <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort")
if (!dir.exists(dest_path_cohort)) {dir.create(dest_path_cohort)}

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 10)

stats <- c("weight_stat", "biomass_stat", "dw_stat", "SGR_stat", "E_somat_stat", 
           "P_excr_stat", "L_excr_stat", "C_excr_stat", "P_uneat_stat", "L_uneat_stat", 
           "C_uneat_stat", "ing_act_stat", "anab_stat", "catab_stat", "NH4_stat", 
           "O2_stat", "food_prov_stat", "rel_feeding_stat", "T_response_stat")
cohort_stats <- c("biomass", "dw", "SGR", "P_excr", "C_excr", "L_excr", "P_uneat", 
                  "C_uneat", "L_uneat", "ing_act", "NH4", "O2", "food_prov")

# Function to process a single farm
process_farm <- function(fid) {
  fnames_farms <- file.path(dest_path_farm, paste0(paste("farmID", fixnum(farm_IDs[fid]), stats, sep = "_"), ".parquet"))
  
  if (overwrite == T | any(!file.exists(fnames_farms))) {
    df_1 <- tar_read(main_farm_growth, branches = fid)
    df_2 <- tar_read(main_farm_growth, branches = length(farm_IDs) + fid)
    df_3 <- tar_read(main_farm_growth, branches = 2*length(farm_IDs) + fid)
    
    walk(seq_along(stats), function(st) {
      if (overwrite == T | !file.exists(fnames_farms[st])) {
        df_11 <- df_1[[st]] %>%
          as.data.frame() %>%
          rename(mean = V1, sd = V2) %>% 
          mutate(feed = feed_types[1], days = 1:548)
        df_21 <- df_2[[st]] %>%
          as.data.frame() %>%
          rename(mean = V1, sd = V2) %>% 
          mutate(feed = feed_types[2], days = 1:548)
        df_31 <- df_3[[st]] %>%
          as.data.frame() %>%
          rename(mean = V1, sd = V2) %>% 
          mutate(feed = feed_types[3], days = 1:548)
        df <- rbind(df_11, df_21, df_31) %>% 
          relocate(days, .before = mean) %>% 
          mutate(farm_ID = farm_IDs[fid], 
                 feed = factor(feed, levels = c("past", "reference", "future"))) %>% 
          write_parquet(fnames_farms[st])
      }
    })
    rm(df_1, df_2, df_3)
  }
  
  fnames_cohorts <- file.path(dest_path_cohort, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts", cohort_stats, sep = "_"), ".parquet"))
  
  if (overwrite == T | any(!file.exists(fnames_cohorts))) {
    df_1 <- tar_read(cohorts_data, branches = fid)
    df_2 <- tar_read(cohorts_data, branches = length(farm_IDs) + fid)
    df_3 <- tar_read(cohorts_data, branches = 2*length(farm_IDs) + fid)
    
    walk(seq_along(cohort_stats), function(st) {
      if (overwrite == T | !file.exists(fnames_cohorts[st])) {
        df <- list(
          df_1[[st]] %>% as.data.frame(),
          df_2[[st]] %>% as.data.frame(),
          df_3[[st]] %>% as.data.frame()
        ) %>% bind_rows() %>% 
          mutate(feed = factor(feed, levels = feed_types)) %>% 
          write_parquet(fnames_cohorts[st])
      }
    })
    rm(df_1, df_2, df_3)
  }
}

# Process all farms in parallel
future_map(seq_along(farm_IDs), process_farm, .progress = TRUE)

# Clean up parallel processing
plan(sequential)
