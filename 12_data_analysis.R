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
  library(furrr)
  library(future)
  conflicts_prefer(dplyr::select(), dplyr::filter(), .quiet = T)
})

# This script aims to take the outputs from the model and answer specific questions 
# Paths & globals -------------------------------------------------------------------------------------------------
source("00_model_functions.R")

coho_path <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort")
proc_path <- file.path("data", "atlantic_salmon", "data_products", "processed_data")

# * Main point is the difference between feeds. 
# * This should be expressed in t/t salmon, so farm biomass is also important.
# * Also interesting to see if difference between feeds varies geographically (correlates with mean/median/max temperature?)

# Biomass produced ------------------------------------------------------------------------------------------------
# How accurately is the model producing the correct biomass production levels?

farm_info <- file.path("data", "_general_data", "farm_locations", "atlantic_salmon_locations_w_temps.qs") %>% 
  qread() %>% 
  filter(day == "day_1") %>% 
  mutate(country = as.factor(country),
         day = str_split_i(day, "_", 2) %>% as.integer()) %>% 
  dplyr::select(-c(model_name, F_CODE, data_year, data_source, details, data_type_2, data_type, species_group, harvest_n, stocking_n, harvest_size_t, daily_mort_rate, day))

hist(farm_info$tonnes_per_farm)

coho_biom_fnms <- file.path(coho_path) %>% 
  list.files(full.names = T) %>% 
  str_subset("biomass")

harv_day <- coho_biom_fnms[1] %>% read_parquet() %>% 
  filter(cohort == 1) %>% 
  distinct(yday) %>% 
  max()
all_days <- coho_biom_fnms[1] %>% read_parquet() %>% 
  distinct(yday) %>% 
  nrow()

coho_biom <- list()
for (i in seq_along(coho_biom_fnms)) {
  df <- coho_biom_fnms[i] %>% 
    read_parquet() %>% 
    filter(yday == harv_day & cohort == 1) %>% 
    select(farm_ID, feed, mean, sd) %>% 
    mutate(sd = sd/mean,
           sd = 3*sd/(all_days/365),
           mean = 3*mean/(all_days/365),
           sd = sd * mean) %>% 
    mutate(mean = mean %>% set_units("g") %>% set_units("t") %>% drop_units(),
           sd = sd %>% set_units("g") %>% set_units("t") %>% drop_units())
  coho_biom[[i]] <- df
}
coho_biom <- bind_rows(coho_biom) %>% 
  merge(farm_info, by.x = "farm_ID", by.y = "farm_id", all.x = T, all.y = F)
qsave(coho_biom, file.path(proc_path, "farm_biomass_produced.qs"))

# Difference between feeds ----------------------------------------------------------------------------------------
plan(multisession, workers = parallel::detectCores() - 10)

farm_IDs <- targets::tar_read(farm_IDs, store = "05_targets_farm")[1:10]

## Total uneaten --------------------------------------------------------------------------------------------------
# Pre-filter files
farm_files <- lapply(farm_IDs, function(id) {
  file.path(coho_path) %>% 
    list.files(full.names = T) %>% 
    str_subset("uneat") %>% 
    str_subset("total", negate = T) %>% 
    str_subset(fixnum(id))
})

# Condense desired processing into a single function
process_uneat <- function(fid, farm_IDs, uneat_fnms) {
  df <- list(
    uneat_fnms %>% str_subset("C_uneat") %>% arrow::read_parquet(mmap = TRUE) %>% mutate(macro = "C"),
    uneat_fnms %>% str_subset("L_uneat") %>% arrow::read_parquet(mmap = TRUE) %>% mutate(macro = "L"),
    uneat_fnms %>% str_subset("P_uneat") %>% arrow::read_parquet(mmap = TRUE) %>% mutate(macro = "P")
  ) %>% 
    bind_rows() %>% 
    filter(prod_day != max(prod_day)) %>% 
    mutate(farm_ID = farm_IDs[fid],
           macro = as.factor(macro),
           stat = "uneaten")
  return(df)
}

# Run with parallel processing
total_uneat <- future_map(
  seq_along(farm_IDs),
  ~process_uneat(.x, farm_IDs, farm_files[[.x]]),
  .progress = T
)
total_uneat <- total_uneat %>% 
  bind_rows() %>% 
  mutate(stat = factor(stat, levels = c("uneaten", "excreted")),
         mean = units::set_units(mean, "g"),
         sd = units::set_units(sd, "g")) %>% 
  arrow::write_parquet(file.path(proc_path, "total_uneaten_raw.parquet"))

# Filter by the end of cohort 1(+1) and the end of cohort 2 - ie where the pattern repeats
ydays <- total_uneat %>% 
  group_by(farm_ID, cohort) %>% 
  reframe(end = max(yday)) %>% 
  pivot_wider(names_from = cohort, names_prefix = "coho_", values_from = end) %>% 
  mutate(start = coho_1+1,
         end = coho_2) %>% 
  select(farm_ID, start, end)

total_uneat <- total_uneat %>% 
  merge(ydays, by = "farm_ID") %>% 
  filter(yday >= start & yday <= end) %>% 
  mutate(sd = sd/mean) %>% 
  group_by(farm_ID, feed, yday, stat, macro) %>% 
  reframe(mean = sum(mean, na.rm = T),
          sd = mean * sum(sd, na.rm = T)) %>% 
  arrow::write_parquet(file.path(proc_path, "total_uneaten.parquet"))

# # Quick plot
# total_uneat %>% 
#   filter(farm_ID == 722) %>% 
#   ggplot(aes(x = yday, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed)) +
#   geom_line(linewidth = 0.75) +
#   geom_ribbon(alpha = 0, linetype = "dotted") +
#   facet_wrap(~macro) +
#   scale_colour_manual(values = c("#666", "#cc0088", "#00aaaa")) +
#   theme_classic()

## Total excreted --------------------------------------------------------------------------------------------------
# Pre-filter files
farm_files <- lapply(farm_IDs, function(id) {
  file.path(coho_path) %>% 
    list.files(full.names = T) %>% 
    str_subset("excr") %>% 
    str_subset("total", negate = T) %>% 
    str_subset(fixnum(id))
})

# Condense desired processing into a single function
process_excr <- function(fid, farm_IDs, excr_fnms) {
  df <- list(
    excr_fnms %>% str_subset("C_excr") %>% arrow::read_parquet(mmap = TRUE) %>% mutate(macro = "C"),
    excr_fnms %>% str_subset("L_excr") %>% arrow::read_parquet(mmap = TRUE) %>% mutate(macro = "L"),
    excr_fnms %>% str_subset("P_excr") %>% arrow::read_parquet(mmap = TRUE) %>% mutate(macro = "P")
  ) %>% 
    bind_rows() %>% 
    filter(prod_day != max(prod_day)) %>% 
    mutate(farm_ID = farm_IDs[fid],
           macro = as.factor(macro),
           stat = "excreted")
  return(df)
}

# Run with parallel processing
total_excr <- future_map(
  seq_along(farm_IDs),
  ~process_excr(.x, farm_IDs, farm_files[[.x]]),
  .progress = T
)
total_excr <- total_excr %>% 
  bind_rows() %>% 
  mutate(stat = factor(stat, levels = c("uneaten", "excreted")),
         mean = units::set_units(mean, "g"),
         sd = units::set_units(sd, "g")) %>% 
  arrow::write_parquet(file.path(proc_path, "total_excreted_raw.parquet"))

# Filter by the end of cohort 1(+1) and the end of cohort 2 - ie where the pattern repeats
ydays <- total_excr %>% 
  group_by(farm_ID, cohort) %>% 
  reframe(end = max(yday)) %>% 
  pivot_wider(names_from = cohort, names_prefix = "coho_", values_from = end) %>% 
  mutate(start = coho_1+1,
         end = coho_2) %>% 
  select(farm_ID, start, end)

total_excr <- total_excr %>% 
  merge(ydays, by = "farm_ID") %>% 
  filter(yday >= start & yday <= end) %>% 
  mutate(sd = sd/mean) %>% 
  group_by(farm_ID, feed, yday, stat, macro) %>% 
  reframe(mean = sum(mean, na.rm = T),
          sd = mean * sum(sd, na.rm = T)) %>% 
  arrow::write_parquet(file.path(proc_path, "total_excreted.parquet"))

plan(sequential)
