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
wt_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "weight_parameter_sensitivity.parquet"))
dw_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "dw_parameter_sensitivity.parquet"))
excr_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "excreted_parameter_sensitivity.parquet"))
uneat_ls %>% bind_rows() %>% write_parquet(file.path("data", "atlantic_salmon", "data_products", "uneaten_parameter_sensitivity.parquet"))

## Farm growth ----------------------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_farm")
# tar_visnetwork()
tar_validate()
# tar_prune()
upd <- tar_outdated(branches = F)
tar_make(reporter = "balanced", seconds_meta_append = 300)

farm_IDs <- tar_read(farm_IDs)
feed_types <- tar_read(feed_types)

overwrite <- T

dest_path_1 <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_farm")
dest_path_2 <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort")

for (fid in 1:length(farm_IDs)) {
  tic()
  ### Process raw outputs -----------------------------------------------------------------------------------------
  # Seperate the different portions of the list - all diets for one farm in each file
  stats <- c("weight_stat", "biomass_stat", "dw_stat", "SGR_stat", "E_somat_stat", "P_excr_stat", "L_excr_stat", "C_excr_stat", "P_uneat_stat", "L_uneat_stat", "C_uneat_stat", "ing_act_stat", "anab_stat", "catab_stat", "NH4_stat", "O2_stat", "food_prov_stat", "rel_feeding_stat", "T_response_stat")
  
  fnames_farms <- file.path(dest_path_1, paste0(paste("farmID", fixnum(farm_IDs[fid]), stats, sep = "_"), ".parquet"))
  
  if (overwrite == T | any(!file.exists(fnames_farms))) {
    df_1 <- tar_read(main_farm_growth, branches = fid)
    df_2 <- tar_read(main_farm_growth, branches = length(farm_IDs) + fid)
    df_3 <- tar_read(main_farm_growth, branches = 2*length(farm_IDs) + fid)
    
    for (st in seq_along(stats)) {
      if (overwrite == T | !file.exists(fnames_farms[st])) {
        df <- list(
          df_1[[st]] %>%
            as.data.frame() %>%
            rename(mean = V1, sd = V2) %>% 
            mutate(feed = feed_types[1], days = 1:548),
          df_2[[st]] %>%
            as.data.frame() %>%
            rename(mean = V1, sd = V2) %>% 
            mutate(feed = feed_types[2], days = 1:548),
          df_3[[st]] %>%
            as.data.frame() %>%
            rename(mean = V1, sd = V2) %>% 
            mutate(feed = feed_types[3], days = 1:548)
        ) %>% bind_rows() %>% 
          relocate(days, .before = mean) %>% 
          mutate(farm_ID = farm_IDs[fid], 
                 feed = factor(feed, levels = c("past", "reference", "future"))) %>% 
          write_parquet(fnames_farms[st])
      }
    }
    rm(df, df_1, df_2, df_3)
  }

  fnames_cohorts <- c(
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_biomass", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_dw", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_SGR", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_P_excr", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_C_excr", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_L_excr", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_P_uneat", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_C_uneat", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_L_uneat", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_ing_act", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_NH4", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_O2", sep = "_"), ".parquet")),
    file.path(dest_path_2, paste0(paste("farmID", fixnum(farm_IDs[fid]), "cohorts_food_prov", sep = "_"), ".parquet"))
  )
  stat_levels <- c("biomass", "dw", "SGR", "P_excr", "C_excr", "L_excr", "P_uneat", "C_uneat", "L_uneat", "ing_act", "NH4", "O2", "food_prov")

  df_1 <- tar_read(cohorts_data, branches = fid)
  df_2 <- tar_read(cohorts_data, branches = length(farm_IDs) + fid)
  df_3 <- tar_read(cohorts_data, branches = 2*length(farm_IDs) + fid)
  
  for (fn in seq_along(fnames_cohorts)) {
    if (overwrite == T | !file.exists(fnames_cohorts[fn])) {
      df <- list(
        df_1[[fn]], 
        df_2[[fn]], 
        df_3[[fn]]
      ) %>% bind_rows() %>% 
        # filter(yday <= 1095 & yday > 365) %>% 
        mutate(feed = factor(feed, levels = feed_types),
               stat = factor(stat, levels = stat_levels)) %>% 
        write_parquet(fnames_cohorts[fn])
    }
  }
  rm(df, df_1, df_2, df_3)
  
  t <- toc(quiet = T)$callback_msg
  print(paste0("Finished processing farm_ID ", farm_IDs[fid], " -- ", round(100*fid/length(farm_IDs),2), "% done -- ", t))
}

# Analysis & plotting ---------------------------------------------------------------------------------------------
## Harvest size ---------------------------------------------------------------------------------------------------
pdata <- tar_read(farm_harvest_size, store = "04_targets_individual")
write_parquet(pdata, file.path("data", "atlantic_salmon", "data_products", "harvest_size.parquet"))

# p <- pdata %>% 
#   mutate(weight = weight %>% set_units("g") %>% set_units("kg")) %>% 
#   select(c(farm_ID, weight)) %>% 
#   ggplot(aes(x = weight)) +
#   geom_histogram(binwidth = 0.25, colour = "black", fill = "salmon", alpha = 0.5) +
#   scale_y_continuous(limits = c(0,500)) +
#   labs(y = "Frequency", x = "Harvest weight") +
#   theme_classic() +
#   theme(legend.position = "none", aspect.ratio = 0.75)
# 
# ggsave(
#   plot = p, 
#   filename = file.path("data/atlantic_salmon/data_products/figures", paste0("harvest_size.png")),
#   height = 150, width = 200, units = "mm"
# )

## Previously processed data --------------------------------------------------------------------------------------
Sys.setenv(TAR_PROJECT = "project_farm")
farm_IDs <- tar_read(farm_IDs)
feed_types <- tar_read(feed_types)
data_path_1 <- file.path(prod_path, "model_outputs_farm")
data_path_2 <- file.path(prod_path, "model_outputs_cohort")
dest_path_1 <- file.path(prod_path, "figures", "farms")
dest_path_2 <- file.path(prod_path, "figures", "cohorts")

make_label <- function(lab){lab %>% str_remove_all("_stat") %>% str_replace_all("_", " ") %>% str_to_title()}

farm_stats <- c("weight_stat", "biomass_stat", "dw_stat", "SGR_stat", "E_somat_stat", "P_excr_stat", "L_excr_stat", "C_excr_stat", "P_uneat_stat", "L_uneat_stat", "C_uneat_stat", "ing_act_stat", "anab_stat", "catab_stat", "NH4_stat", "O2_stat", "food_prov_stat", "rel_feeding_stat", "T_response_stat")

# List all the files output from targets processing
farm_files <- list()
for (s in seq_along(farm_stats)) {
  farm_files[[s]] <- list.files(data_path_1, full.names = T) %>% 
    str_subset(farm_stats[s])
}

## Save farm seperate plots ---------------------------------------------------------------------------------------
# for (ff in seq_along(farm_files)) {
#   fig_path <- file.path(dest_path_1, farm_stats[[ff]])
#   if (!dir.exists(fig_path)) {dir.create(fig_path)}
#   
#   for (fid in seq_along(farm_IDs)) {
#     png_name <- file.path(fig_path, paste0(farm_stats[ff], "_", fixnum(farm_IDs[fid]), ".png"))
#     
#     if (overwrite == T | !file.exists(png_name)) {
#       df <- farm_files[[ff]] %>%
#         str_subset(fixnum(farm_IDs[fid])) %>%
#         read_parquet()
#       
#       if (farm_stats[ff] %in% c("weight_stat", "biomass_stat")) {
#         df$mean = set_units(set_units(df$mean, "g"), "kg")
#         df$sd   = set_units(set_units(df$sd, "g"), "kg")
#       } else if (farm_stats[ff] %in% c("E_somat_stat")) {
#         df$mean = set_units(df$mean, "J")
#         df$sd   = set_units(df$sd, "J")
#       } else if (farm_stats[ff] %in% c("anab_stat", "catab_stat")) {
#         df$mean = set_units(df$mean, "J d-1")
#         df$sd   = set_units(df$sd, "J d-1")
#       } else if (farm_stats[ff] %in% c("SGR_stat", "rel_feeding_stat", "T_response_stat")) {
#         df$mean = df$mean
#         df$sd   = df$sd
#       } else {
#         df$mean = set_units(df$mean, "g d-1")
#         df$sd   = set_units(df$sd, "g d-1")
#       }
#       
#       p <- df %>% 
#         ggplot(aes(x = days, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed, fill = feed)) +
#         geom_line(linewidth = 0.75) +
#         geom_ribbon(alpha = 0, linetype = "dashed") +
#         scale_colour_brewer(palette = "Set1") +
#         scale_fill_brewer(palette = "Set1") +
#         scale_x_continuous(breaks = seq(0, 550, 100)) +
#         theme_classic() +
#         theme(aspect.ratio = 0.75, legend.position = "none") +
#         labs(x = "Day of production", y = make_label(farm_stats[ff]))
#         
#       ggsave(
#         plot = p, 
#         filename = png_name,
#         height = 150, width = 200, units = "mm"
#       )
#     }
#   }
#   print(paste0("Finished processing all farms, stat ", ff, " -- ", round(100*ff/length(farm_files),2), "% done"))
# }

## Combine some stats ---------------------------------------------------------------------------------------------
### Farms ---------------------------------------------------------------------------------------------------------
farm_files <- farm_files %>% 
  Reduce(c, .) %>% 
  str_subset(".parquet")

# Combined excretion
print("Starting farm total excretion processing...")
path_2 <- file.path("data", "atlantic_salmon", "data_products", "figures", "farms", "total_excr_stat")
if (!dir.exists(path_2)) {dir.create(path_2)}
for (fid in seq_along(farm_IDs)) {
  # Total excretion (carbs, lipis, and fats)
  par_name <- file.path(data_path_1, paste0("total_excr_stat_", fixnum(farm_IDs[fid]), ".parquet"))
  png_name <- file.path(path_2, paste0("total_excr_stat_", fixnum(farm_IDs[fid]), ".png"))
  
  if (overwrite == T | any(!file.exists(par_name), !file.exists(png_name))) {
    fnms <- farm_files %>% 
      str_subset(fixnum(farm_IDs[fid])) %>% 
      str_subset("excr_stat")
    
    df <- list(
      fnms %>% str_subset("C_excr") %>% read_parquet() %>% mutate(stat = "C_excr"),
      fnms %>% str_subset("L_excr") %>% read_parquet() %>% mutate(stat = "L_excr"),
      fnms %>% str_subset("P_excr") %>% read_parquet() %>% mutate(stat = "P_excr")
    ) %>% 
      bind_rows() %>% 
      filter(days != max(days)) %>% 
      write_parquet(par_name)
    
    # p <- df %>% 
    #   group_by(days, feed, farm_ID) %>% 
    #   reframe(mean = sum(mean, na.rm = T),
    #           sd = sum(sd, na.rm = T)) %>% 
    #   mutate(mean = set_units(mean, "g d-1") %>% set_units("kg d-1"),
    #          sd = set_units(sd, "g d-1") %>% set_units("kg d-1")) %>% 
    #   ggplot(aes(x = days, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed, fill = feed)) +
    #   geom_line(linewidth = 0.75) +
    #   geom_ribbon(alpha = 0, linetype = "dashed") +
    #   scale_colour_brewer(palette = "Set1") +
    #   scale_fill_brewer(palette = "Set1") +
    #   scale_x_continuous(breaks = seq(0, 550, 100)) +
    #   theme_classic() +
    #   theme(aspect.ratio = 0.75, legend.position = "none") +
    #   labs(x = "Day of production", y = "Total excretion")
    # 
    # ggsave(
    #   plot = p, 
    #   filename = png_name,
    #   height = 150, width = 200, units = "mm"
    # )
  }
  if (fid %in% as.integer(seq(0, length(farm_IDs), length.out = 10))) {
    print(paste0("Finish processing ", fid, " of ", length(farm_IDs), " -- ", round(100*fid/length(farm_IDs),2), "% done"))
  }
}

# Combined uneaten
print("Starting farm total uneaten feed processing...")
path_2 <- file.path("data", "atlantic_salmon", "data_products", "figures", "farms", "total_uneat_stat")
if (!dir.exists(path_2)) {dir.create(path_2)}
for (fid in seq_along(farm_IDs)) {
  # Total uneaten (carbs, lipis, and fats)
  par_name <- file.path(data_path_1, paste0("total_excr_stat_", fixnum(farm_IDs[fid]), ".parquet"))
  png_name <- file.path(path_2, paste0("total_excr_stat_", fixnum(farm_IDs[fid]), ".png"))
  
  if (overwrite == T | any(!file.exists(par_name), !file.exists(png_name))) {
    fnms <- farm_files %>% 
      str_subset(fixnum(farm_IDs[fid])) %>% 
      str_subset("uneat_stat")
    
    df <- list(
      fnms %>% str_subset("C_uneat") %>% read_parquet() %>% mutate(stat = "C_uneat"),
      fnms %>% str_subset("L_uneat") %>% read_parquet() %>% mutate(stat = "L_uneat"),
      fnms %>% str_subset("P_uneat") %>% read_parquet() %>% mutate(stat = "P_uneat")
    ) %>% 
      bind_rows() %>% 
      filter(days != max(days)) %>% 
      mutate(stat = as.factor(stat)) %>% 
      write_parquet(par_name)
    
    # p <- df %>% 
    #   group_by(days, feed, farm_ID) %>% 
    #   reframe(mean = sum(mean, na.rm = T),
    #           sd = sum(sd, na.rm = T)) %>% 
    #   mutate(mean = set_units(mean, "g d-1") %>% set_units("kg d-1"),
    #          sd = set_units(sd, "g d-1") %>% set_units("kg d-1")) %>% 
    #   ggplot(aes(x = days, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed, fill = feed)) +
    #   geom_line(linewidth = 0.75) +
    #   geom_ribbon(alpha = 0, linetype = "dashed") +
    #   scale_colour_brewer(palette = "Set1") +
    #   scale_fill_brewer(palette = "Set1") +
    #   scale_x_continuous(breaks = seq(0, 550, 100)) +
    #   theme_classic() +
    #   theme(aspect.ratio = 0.75, legend.position = "none") +
    #   labs(x = "Day of production", y = "Total uneaten")
    # 
    # ggsave(
    #   plot = p, 
    #   filename = png_name,
    #   height = 150, width = 200, units = "mm"
    # )
  }
  if (fid %in% as.integer(seq(0, length(farm_IDs), length.out = 10))) {
    print(paste0("Finish processing ", fid, " of ", length(farm_IDs), " -- ", round(100*fid/length(farm_IDs),2), "% done"))
  }
}

### Cohorts -------------------------------------------------------------------------------------------------------
cohort_stats <- c("cohorts_biomass", "cohorts_dw", "cohorts_SGR", "cohorts_O2", "cohorts_P_uneat", "cohorts_C_uneat", "cohorts_L_uneat", "cohorts_food_prov", "cohorts_NH4", "cohorts_P_excr", "cohorts_C_excr", "cohorts_L_excr")

cohort_files <- list()
for (s in seq_along(cohort_stats)) {
  cohort_files[[s]] <- list.files(data_path_2, full.names = T) %>% str_subset(cohort_stats[s])
}
cohort_files <- cohort_files %>% 
  Reduce(c, .) %>% 
  str_subset(".parquet")

# Combined excretion
print("Starting cohort total excretion processing...")
path_2 <- file.path("data", "atlantic_salmon", "data_products", "figures", "cohorts", "total_excr_stat")
if (!dir.exists(path_2)) {dir.create(path_2)}
for (fid in seq_along(farm_IDs)) {
  # Total excretion (carbs, lipis, and fats)
  par_name <- file.path(data_path_2, paste0("total_excr_stat_", fixnum(farm_IDs[fid]), ".parquet"))
  png_name <- file.path(path_2, paste0("total_excr_stat_", fixnum(farm_IDs[fid]), ".png"))
  
  if (overwrite == T | any(!file.exists(par_name), !file.exists(png_name))) {
    fnms <- cohort_files %>% 
      str_subset(fixnum(farm_IDs[fid])) %>% 
      str_subset("excr")
    
    df <- list(
      fnms %>% str_subset("C_excr") %>% read_parquet() %>% mutate(stat = "C_excr"),
      fnms %>% str_subset("L_excr") %>% read_parquet() %>% mutate(stat = "L_excr"),
      fnms %>% str_subset("P_excr") %>% read_parquet() %>% mutate(stat = "P_excr")
    ) %>% 
      bind_rows() %>% 
      filter(yday != max(yday)) %>% 
      mutate(stat = as.factor(stat),
             cohort = as.factor(cohort)) %>% 
      write_parquet(par_name)
    
    # p <- df %>% 
    #   group_by(yday, feed, farm_ID) %>% 
    #   reframe(mean = sum(mean, na.rm = T),
    #           sd = sum(sd, na.rm = T)) %>% 
    #   mutate(mean = set_units(mean, "g d-1") %>% set_units("kg d-1"),
    #          sd = set_units(sd, "g d-1") %>% set_units("kg d-1")) %>% 
    #   ggplot(aes(x = yday, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed, fill = feed)) +
    #   geom_line(linewidth = 0.75) +
    #   geom_ribbon(alpha = 0, linetype = "dashed") +
    #   scale_colour_brewer(palette = "Set1") +
    #   scale_fill_brewer(palette = "Set1") +
    #   scale_x_continuous(breaks = seq(0, 1600, 200)) +
    #   theme_classic() +
    #   theme(aspect.ratio = 0.75, legend.position = "none") +
    #   labs(x = "Day of production", y = "Total excretion")
    # 
    # ggsave(
    #   plot = p, 
    #   filename = png_name,
    #   height = 150, width = 200, units = "mm"
    # )
  }
  if (fid %in% as.integer(seq(0, length(farm_IDs), length.out = 10))) {
    print(paste0("Finish processing ", fid, " of ", length(farm_IDs), " -- ", round(100*fid/length(farm_IDs),2), "% done"))
  }
}

# Combined uneaten
print("Starting cohort total uneaten feed processing...")
path_2 <- file.path("data", "atlantic_salmon", "data_products", "figures", "cohorts", "total_uneat_stat")
if (!dir.exists(path_2)) {dir.create(path_2)}
for (fid in seq_along(farm_IDs)) {
  # Total uneaten (carbs, lipis, and fats)
  par_name <- file.path(data_path_2, paste0("total_excr_stat_", fixnum(farm_IDs[fid]), ".parquet"))
  png_name <- file.path(path_2, paste0("total_excr_stat_", fixnum(farm_IDs[fid]), ".png"))
  
  if (overwrite == T | any(!file.exists(par_name), !file.exists(png_name))) {
    fnms <- cohort_files %>% 
      str_subset(fixnum(farm_IDs[fid])) %>% 
      str_subset("uneat")
    
    df <- list(
      fnms %>% str_subset("C_uneat") %>% read_parquet() %>% mutate(stat = "C_uneat"),
      fnms %>% str_subset("L_uneat") %>% read_parquet() %>% mutate(stat = "L_uneat"),
      fnms %>% str_subset("P_uneat") %>% read_parquet() %>% mutate(stat = "P_uneat")
    ) %>% 
      bind_rows() %>% 
      filter(yday != max(yday)) %>% 
      mutate(stat = as.factor(stat),
             cohort = as.factor(cohort)) %>% 
      write_parquet(par_name)
    
    # p <- df %>% 
    #   group_by(yday, feed, farm_ID) %>% 
    #   reframe(mean = sum(mean, na.rm = T),
    #           sd = sum(sd, na.rm = T)) %>% 
    #   mutate(mean = set_units(mean, "g d-1") %>% set_units("kg d-1"),
    #          sd = set_units(sd, "g d-1") %>% set_units("kg d-1")) %>% 
    #   ggplot(aes(x = yday, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed, fill = feed)) +
    #   geom_line(linewidth = 0.75) +
    #   geom_ribbon(alpha = 0, linetype = "dashed") +
    #   scale_colour_brewer(palette = "Set1") +
    #   scale_fill_brewer(palette = "Set1") +
    #   scale_x_continuous(breaks = seq(0, 550, 100)) +
    #   theme_classic() +
    #   theme(aspect.ratio = 0.75, legend.position = "none") +
    #   labs(x = "Day of production", y = "Total uneaten")
    # 
    # ggsave(
    #   plot = p, 
    #   filename = png_name,
    #   height = 150, width = 200, units = "mm"
    # )
  }
  if (fid %in% as.integer(seq(0, length(farm_IDs), length.out = 10))) {
    print(paste0("Finish processing ", fid, " of ", length(farm_IDs), " -- ", round(100*fid/length(farm_IDs),2), "% done"))
  }
}

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
