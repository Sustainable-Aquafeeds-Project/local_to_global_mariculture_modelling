# Setup -----------------------------------------------------------------------------------------------------------
# library(tidyverse)
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
library(conflicted)
conflicts_prefer(dplyr::select(), dplyr::filter(), .quiet = T)

# functions
source("00_model_functions.R")

# species paths
this_species <- "atlantic_salmon"
this_path <- sprintf("data/%s", this_species)

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
  filter(mean_temp < 4) %>% 
  pull(farm_id)

qsave(x = farms_to_omit, 
      file = sprintf("data/_general_data/farm_locations/%s_farms_to_omit.qs", this_species))


# STEP 2 - Run model ----------------------------------------------------------------------------------------------
## Example individuals --------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_individual")
# tar_visnetwork()
tar_make(reporter = "summary", seconds_meta_append = 300)
tar_prune()

farm_IDs <- tar_read(farm_IDs)
sens_all_params <- tar_read(sens_all_params)
patt <- data.frame(
  farm_ID = rep(farm_IDs, each = length(sens_all_params)*3),
  param = rep(rep(1:length(sens_all_params), each = 3), times = length(farm_IDs)),
  factor = rep(c(0.9,1,1.1), times = length(sens_all_params)*length(farm_IDs))
)
patt$br <- 1:nrow(patt)

wt_ls <- dw_ls <- excr_ls <- uneat_ls <- list()
for (p in 1:length(sens_all_params)) {
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
wt_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "weight_parameter_sensitivity.parquet"))
dw_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "dw_parameter_sensitivity.parquet"))
excr_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "excreted_parameter_sensitivity.parquet"))
uneat_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "uneaten_parameter_sensitivity.parquet"))

## Farm growth ----------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_farm")
# tar_visnetwork()
tar_make(reporter = "summary", seconds_meta_append = 300)
# tar_prune()

farm_IDs <- tar_read(farm_IDs)
feed_types <- tar_read(feed_types)
stats <- c("weight", "biomass", "dw", "SGR", "E_somat", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "ing_act", "anab", "catab", "O2", "NH4", "food_prov", "rel_feeding", "T_response", "total_excr_mat", "total_uneat_mat")

overwrite <- F

for (f in 1:length(farm_IDs)) {
  for (i in 1:length(stats)) {
    fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_farm", paste0(paste("farmID", fixnum(farm_IDs[f]), stats[i], sep = "_"), ".parquet"))
    if (overwrite == T | !file.exists(fname)) {
      df <- list(
        tar_read(main_farm_growth, branches = f)[[i]] %>%
          as.data.frame() %>%
          rename(mean = V1, sd = V2) %>% 
          mutate(feed = feed_types[1], days = 1:548),
        tar_read(main_farm_growth, branches = length(farm_IDs) + f)[[i]] %>%
          as.data.frame() %>%
          rename(mean = V1, sd = V2) %>% 
          mutate(feed = feed_types[2], days = 1:548),
        tar_read(main_farm_growth, branches = 2*length(farm_IDs) + f)[[i]] %>%
          as.data.frame() %>%
          rename(mean = V1, sd = V2) %>% 
          mutate(feed = feed_types[3], days = 1:548)
      ) %>% bind_rows() %>% 
        relocate(days, .before = mean) %>% 
        mutate(farm_ID = farm_IDs[f], feed = factor(feed, levels = feed_types)) %>% 
        write_parquet(fname)
    }
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_biomass", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_biomass, branches = f), 
      tar_read(cohorts_biomass, branches = length(farm_IDs) + f), 
      tar_read(cohorts_biomass, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
     write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_dw", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_dw, branches = f), 
      tar_read(cohorts_dw, branches = length(farm_IDs) + f), 
      tar_read(cohorts_dw, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_SGR", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_SGR, branches = f), 
      tar_read(cohorts_SGR, branches = length(farm_IDs) + f), 
      tar_read(cohorts_SGR, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_O2", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_O2, branches = f), 
      tar_read(cohorts_O2, branches = length(farm_IDs) + f), 
      tar_read(cohorts_O2, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }

  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_P_uneat", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_P_uneat, branches = f), 
      tar_read(cohorts_P_uneat, branches = length(farm_IDs) + f), 
      tar_read(cohorts_P_uneat, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_C_uneat", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_C_uneat, branches = f), 
      tar_read(cohorts_C_uneat, branches = length(farm_IDs) + f), 
      tar_read(cohorts_C_uneat, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_total_uneat", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_total_uneat_mat, branches = f), 
      tar_read(cohorts_total_uneat_mat, branches = length(farm_IDs) + f), 
      tar_read(cohorts_total_uneat_mat, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_ing_act", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_ing_act, branches = f), 
      tar_read(cohorts_ing_act, branches = length(farm_IDs) + f), 
      tar_read(cohorts_ing_act, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_food_prov", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_food_prov, branches = f), 
      tar_read(cohorts_food_prov, branches = length(farm_IDs) + f), 
      tar_read(cohorts_food_prov, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_NH4", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_NH4, branches = f), 
      tar_read(cohorts_NH4, branches = length(farm_IDs) + f), 
      tar_read(cohorts_NH4, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_L_uneat", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_L_uneat, branches = f), 
      tar_read(cohorts_L_uneat, branches = length(farm_IDs) + f), 
      tar_read(cohorts_L_uneat, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_P_excr", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_P_excr, branches = f), 
      tar_read(cohorts_P_excr, branches = length(farm_IDs) + f), 
      tar_read(cohorts_P_excr, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_L_excr", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_L_excr, branches = f), 
      tar_read(cohorts_L_excr, branches = length(farm_IDs) + f), 
      tar_read(cohorts_L_excr, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_C_excr", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_C_excr, branches = f), 
      tar_read(cohorts_C_excr, branches = length(farm_IDs) + f), 
      tar_read(cohorts_C_excr, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  fname <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort", paste0(paste("farmID", fixnum(farm_IDs[f]), "cohorts_total_excr", sep = "_"), ".parquet"))
  if (overwrite == T | !file.exists(fname)) {
    df <- list(
      tar_read(cohorts_total_excr_mat, branches = f), 
      tar_read(cohorts_total_excr_mat, branches = length(farm_IDs) + f), 
      tar_read(cohorts_total_excr_mat, branches = 2*length(farm_IDs) + f)
    ) %>% bind_rows() %>% 
      filter(days <= 1095 & days > 365) %>% 
      mutate(feed = factor(feed, levels = feed_types)) %>% 
      write_parquet(fname)
  }
  
  print(paste("Farm", f, "of", length(farm_IDs), "(", farm_IDs[f], ")", " finished and saved at", Sys.time()))
}

# Analysis & plotting ---------------------------------------------------------------------------------------------
data_path <- file.path("data", "atlantic_salmon", "data_products")

## Farm production ------------------------------------------------------------------------------------------------
biom_farms <- file.path(data_path, "model_outputs_farm") %>% 
  str_subset("biomass")

## Feed diffs -----------------------------------------------------------------------------------------------------

## Plots ----------------------------------------------------------------------------------------------------------
### Example individuals -------------------------------------------------------------------------------------------
files <- file.path("data_processed", "example_individual", "raw") %>% list.files(full.names = T)
p0 <- read_parquet(files[1]) %>% 
  ggplot(aes(x = days, y = value, colour = feed)) +
  geom_line(linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, 730, by = 182)) +
  scale_colour_brewer(palette = "Set2") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "sans", colour = "black", size = 12))

overwrite <- T

# Not plotted:
# food_enc
# ing_pot
# rel_feeding
# E_somat

max <- tar_read(ind_max_weight) %>% unname() %>% set_units("g") %>% set_units("kg") %>% drop_units() %>% ceiling()
min <- tar_read(ind_min_dw) %>% unname() %>% set_units("g") %>% set_units("kg") %>% drop_units() %>% floor()
weight_fnms <- str_subset(files, "weight")
for (i in 1:length(weight_fnms)){
  fname <- weight_fnms[i] %>% str_replace_all("raw", "plots/weight") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(weight_fnms[i])
    df$value <- set_units(df$value, "g") %>% set_units("kg") %>% drop_units()
    p <- p0 %+% df + 
      scale_y_continuous(breaks = seq(0, 8, 2), limits = c(min, max)) +
      labs(x = "Day of production", y = "Individual weight (kg)")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  
    print(paste(i, "of", length(weight_fnms), "plots saved"))
  }
}

max <- tar_read(ind_max_dw) %>% unname() %>% ceiling()
min <- tar_read(ind_min_dw) %>% unname() %>% floor()
dw_fnms <- str_subset(files, "dw")
for (i in 1:length(dw_fnms)){
  fname <- dw_fnms[i] %>% str_replace_all("raw", "plots/dw") %>% str_replace_all(".parquet", ".png")
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(dw_fnms[i])
    p <- p0 %+% df + 
      scale_y_continuous(breaks = seq(-6, 12, 2), limits = c(-6,max)) +
      geom_hline(aes(yintercept = 0), linetype = "dotted") +
      labs(x = "Day of production", y = expression("Weight change (g d"^-1*")"))
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
    
    print(paste(i, "of", length(dw_fnms), "plots saved"))
  }
}

max <- tar_read(farm_max_T_response) %>% unname() %>% max() %>% set_units("g") %>% set_units("kg") %>% drop_units() %>% ceiling()
Tresp_fnms <- str_subset(files, "T_response")
for (i in 1:length(Tresp_fnms)){
  fname <- Tresp_fnms[i] %>% str_replace_all("raw", "plots/T_response") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(Tresp_fnms[i])
    p <- p0 %+% df + 
      scale_y_continuous(breaks = seq(0, 4, 0.5), limits = c(1, 2.75)) +
      labs(x = "Day of production", y = "Temperature response")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

Pexcr_fnms <- str_subset(files, "P_excr")
for (i in 1:length(Pexcr_fnms)){
  fname <- Pexcr_fnms[i] %>% str_replace_all("raw", "plots/P_excr") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(Pexcr_fnms[i])
    p <- p0 %+% df + 
      scale_y_continuous(breaks = seq(0, 4, 0.2), limits = c(0, 1)) +
      labs(x = "Day of production", y = expression("Protein excreted (g d"^-1*")"))
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

Lexcr_fnms <- str_subset(files, "L_excr")
for (i in 1:length(Lexcr_fnms)){
  fname <- Lexcr_fnms[i] %>% str_replace_all("raw", "plots/L_excr") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(Lexcr_fnms[i])
    p <- p0 %+% df + 
      scale_y_continuous(breaks = seq(0, 4, 0.05), limits = c(0, 0.2)) +
      labs(x = "Day of production", y = expression("Lipid excreted (g d"^-1*")"))
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

Cexcr_fnms <- str_subset(files, "C_excr")
for (i in 1:length(Cexcr_fnms)){
  fname <- Cexcr_fnms[i] %>% str_replace_all("raw", "plots/C_excr") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(Cexcr_fnms[i])
    p <- p0 %+% df + 
      scale_y_continuous(breaks = seq(0, 4, 0.2), limits = c(0, 1.2)) +
      labs(x = "Day of production", y = expression("Carbohydrate excreted (g d"^-1*")"))
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

Puneat_fnms <- str_subset(files, "P_uneat")
for (i in 1:length(Puneat_fnms)){
  fname <- Puneat_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(Puneat_fnms[i])
    df$value <- set_units(df$value, "g d-1") %>% drop_units()
    p <- p0 %+% df + labs(x = "Day of production", y = "Protein uneaten")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

Luneat_fnms <- str_subset(files, "L_uneat")
for (i in 1:length(Luneat_fnms)){
  fname <- Luneat_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(Luneat_fnms[i])
    df$value <- set_units(df$value, "g d-1") %>% drop_units()
    p <- p0 %+% df + labs(x = "Day of production", y = "Lipid uneaten")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

Cuneat_fnms <- str_subset(files, "C_uneat")
for (i in 1:length(Cuneat_fnms)){
  fname <- Cuneat_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(Cuneat_fnms[i])
    df$value <- set_units(df$value, "g d-1") %>% drop_units()
    p <- p0 %+% df + labs(x = "Day of production", y = "Carbohydrate uneaten")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

prov_fnms <- str_subset(files, "food_prov")
for (i in 1:length(prov_fnms)){
  fname <- prov_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(prov_fnms[i])
    df$value <- set_units(df$value, "g d-1") %>% drop_units()
    p <- p0 %+% df + labs(x = "Day of production", y = "Food provided")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

act_fnms <- str_subset(files, "ing_act")
for (i in 1:length(act_fnms)){
  fname <- act_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(act_fnms[i])
    df$value <- set_units(df$value, "g d-1") %>% drop_units()
    p <- p0 %+% df + labs(x = "Day of production", y = "Food ingested")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

assim_fnms <- str_subset(files, "E_assim")
for (i in 1:length(assim_fnms)){
  fname <- assim_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(assim_fnms[i])
    df$value <- set_units(df$value, "J d-1") %>% drop_units()
    p <- p0 %+% df + labs(x = "Day of production", y = "Energy assimilated")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

O2_fnms <- str_subset(files, "O2")
for (i in 1:length(O2_fnms)){
  fname <- O2_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(O2_fnms[i])
    df$value <- set_units(df$value, "g d-1") %>% drop_units()
    p <- p0 %+% df + labs(x = "Day of production", y = "Oxygen consumption")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

NH4_fnms <- str_subset(files, "NH4")
for (i in 1:length(NH4_fnms)){
  fname <- NH4_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(NH4_fnms[i])
    df$value <- set_units(df$value, "g d-1") %>% drop_units()
    p <- p0 %+% df + labs(x = "Day of production", y = expression("NH"_4*"produced"))
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

SGR_fnms <- str_subset(files, "SGR")
for (i in 1:length(SGR_fnms)){
  fname <- SGR_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(SGR_fnms[i])
    p <- p0 %+% df + labs(x = "Day of production", y = "SGR")
    ggsave(plot = p, filename = fname, height = 4.5, width = 7)
  }
}

### Farms ---------------------------------------------------------------------------------------------------------
files <- file.path("data_processed", "farm_growth", "raw") %>% list.files(full.names = T)
p0 <- read_parquet(files[1]) %>% 
  ggplot(aes(x = days, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed, fill = feed)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(alpha = 0.25) + 
  theme_classic() +
  theme(text = element_text(family = "sans", colour = "black", size = 12))

weight_fnms <- str_subset(files, "weight")
for (i in 1:length(weight_fnms)){
  fname <- weight_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(weight_fnms[i])
    df$mean <- set_units(df$mean, "g") %>% set_units("kg")
    df$sd <- set_units(df$sd, "g") %>% set_units("kg")
    p <- p0 %+% df + labs(x = "Day of production", y = "Individual weight")
    ggsave(plot = p, filename = fname)
  }
}

biom_fnms <- str_subset(files, "biomass")
for (i in 1:length(biom_fnms)){
  fname <- biom_fnms[i] %>% str_replace_all("raw", "plots") %>% str_replace_all(".parquet", ".png")
  
  if (overwrite == T | !file.exists(fname)) {
    df <- read_parquet(biom_fnms[i])
    df$mean <- set_units(df$mean, "g") %>% set_units("t")
    df$sd <- set_units(df$sd, "g") %>% set_units("t")
    p <- p0 %+% df + labs(x = "Day of production", y = "Farm biomass")
    ggsave(plot = p, filename = fname)
  }
}



# Old stuff -------------------------------------------------------------------------------------------------------
# STEP 2b - Run model for each location under past feed
feed_type <- "past"
farm_list <- 
  qread(file = sprintf("data/_general_data/farm_locations/%s_locations_w_temps.qs", this_species)) %>% 
  filter(model_name == this_species) %>% 
  group_by(farm_id) %>% 
  group_split()

future::plan(strategy = "multisession", workers = parallel::detectCores()-2) #select core cluster

furrr::future_map(.x = farm_list, .f = \(this_farm){
  this_stocking_N = unique(round(this_farm$stocking_n))
  this_farm_id <- unique(paste0("farmID_", this_farm$farm_id))
  out_loader <- data_loader(Path =  this_path, Farm_id = this_farm_id)
  model_run(Path = this_path, Forcings = out_loader, Feed_type = feed_type, Stocking = this_stocking_N, Farm_id = this_farm_id)
},
.options = furrr_options(seed = 123))

# STEP 2c - Run model for each location under future feed
feed_type <- "future"
farm_list <- 
  qread(file = sprintf("data/_general_data/farm_locations/%s_locations_w_temps.qs", this_species)) %>% 
  filter(model_name == this_species) %>% 
  group_by(farm_id) %>% 
  group_split()

future::plan(strategy = "multisession", workers = parallel::detectCores()-2) #select core cluster

furrr::future_map(.x = farm_list, .f = \(this_farm){
  this_stocking_N = unique(round(this_farm$stocking_n))
  this_farm_id <- unique(paste0("farmID_", this_farm$farm_id))
  out_loader <- data_loader(Path =  this_path, Farm_id = this_farm_id)
  model_run(Path = this_path, Forcings = out_loader, Feed_type = feed_type, Stocking = this_stocking_N, Farm_id = this_farm_id)
},
.options = furrr_options(seed = 123))

# Benchmarking ----------------------------------------------------------------------------------------------------
rbenchmark::benchmark(
  df = {apportion_feed(provided = 10, ingested = 9,
                           prop = feed_params[['Carbohydrates']]$proportion,
                           macro = feed_params[['Carbohydrates']]$ing_carb,
                           digestibility = feed_params[['Carbohydrates']]$ing_carb_digestibility)},
  matrix = {apportion_feed_short(provided = 10, ingested = 9,
                             prop = feed_params[['Carbohydrates']]$proportion,
                             macro = feed_params[['Carbohydrates']]$ing_carb,
                             digestibility = feed_params[['Carbohydrates']]$ing_carb_digestibility)},
  replications = 1000
)
