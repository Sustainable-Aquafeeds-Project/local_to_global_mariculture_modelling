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
if(!dir.exists(proc_path)) {dir.create(proc_path)}

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

coho_biom <- list()
for (i in seq_along(coho_biom_fnms)) {
  df <- coho_biom_fnms[i] %>% 
    read_parquet() %>% 
    select(yday, farm_ID, feed, mean, sd, cohort) %>% 
    mutate(sd = sd/mean,
           sd = 3*sd/(all_days/365),
           mean = 3*mean/(all_days/365),
           sd = sd * mean) %>% 
    mutate(mean = mean %>% set_units("g") %>% set_units("t") %>% drop_units(),
           sd = sd %>% set_units("g") %>% set_units("t") %>% drop_units())
  coho_biom[[i]] <- df
}
coho_biom <- bind_rows(coho_biom)

ydays <- coho_biom %>% 
  group_by(farm_ID, cohort) %>% 
  reframe(end = max(yday)) %>% 
  pivot_wider(names_from = cohort, names_prefix = "coho_", values_from = end) %>% 
  mutate(start = coho_1+1,
         end = coho_2) %>% 
  select(farm_ID, start, end)

coho_biom <- coho_biom %>% 
  merge(ydays, by = "farm_ID") %>% 
  dplyr::filter(yday >= start & yday < end) %>% 
  merge(farm_info, by.x = "farm_ID", by.y = "farm_id", all.x = T, all.y = F)
qsave(coho_biom, file.path(proc_path, "farm_biomass_produced.qs"))

# Difference between feeds ----------------------------------------------------------------------------------------
## Total uneaten --------------------------------------------------------------------------------------------------
uneat_fnms <- list.files("data/atlantic_salmon/data_products/processed_cohort", full.names = T) %>% str_subset("total_uneat")
uneat <- list()
for (f in seq_along(excr_fnms)) {
  uneat[[f]] <- uneat_fnms[f] %>% arrow::read_parquet()
}
uneat <- uneat %>% 
  bind_rows() %>% 
  mutate(macro = str_split_i(stat, "_", 1),
         stat = str_split_i(stat, "_", 2)) %>% 
  mutate(macro = as.factor(macro),
         stat = factor(stat, levels = c("uneat", "excr")))

# Filter by the end of cohort 1(+1) and the end of cohort 2 - ie where the pattern repeats
# Begins after the harvesting of one cohort (and the second cohort is already underway) and ends at the harvesting of the second cohort (when the third cohort is already out)
ydays <- uneat %>% 
  group_by(farm_ID, cohort) %>% 
  reframe(end = max(yday)) %>% 
  pivot_wider(names_from = cohort, names_prefix = "coho_", values_from = end) %>% 
  mutate(start = coho_1+1,
         end = coho_2) %>% 
  select(farm_ID, start, end)

uneat <- uneat %>% 
  merge(ydays, by = "farm_ID") %>% 
  dplyr::filter(yday >= start & yday < end) %>% 
  mutate(sd = sd/mean) %>% 
  group_by(farm_ID, feed, yday, stat, macro) %>% 
  reframe(mean = sum(mean, na.rm = T),
          sd = mean * sum(sd, na.rm = T)) %>% 
  arrow::write_parquet(file.path(proc_path, "total_uneat_allfarms.parquet"))

# # Quick plot
# uneat %>%
#   dplyr::filter(farm_ID == 722) %>%
#   ggplot(aes(x = yday, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed)) +
#   geom_line(linewidth = 0.75) +
#   geom_ribbon(alpha = 0, linetype = "dotted") +
#   facet_wrap(~macro) +
#   # scale_colour_manual(values = c("#666", "#cc0088", "#00aaaa")) +
#   theme_classic()

rm(uneat)

## Total excreted --------------------------------------------------------------------------------------------------
excr_fnms <- list.files("data/atlantic_salmon/data_products/processed_cohort", full.names = T) %>% str_subset("total_excr")
excr <- list()
for (f in seq_along(excr_fnms)) {
  excr[[f]] <- excr_fnms[f] %>% arrow::read_parquet()
}
excr <- excr %>% 
  bind_rows() %>% 
  mutate(macro = str_split_i(stat, "_", 1),
         stat = str_split_i(stat, "_", 2)) %>% 
  mutate(macro = as.factor(macro),
         stat = factor(stat, levels = c("uneat", "excr")))

# Filter by the end of cohort 1(+1) and the end of cohort 2 - ie where the pattern repeats
ydays <- excr %>% 
  group_by(farm_ID, cohort) %>% 
  reframe(end = max(yday)) %>% 
  pivot_wider(names_from = cohort, names_prefix = "coho_", values_from = end) %>% 
  mutate(start = coho_1+1,
         end = coho_2) %>% 
  select(farm_ID, start, end)

excr <- excr %>% 
  merge(ydays, by = "farm_ID") %>% 
  dplyr::filter(yday >= start & yday < end) %>% 
  mutate(sd = sd/mean) %>% 
  group_by(farm_ID, feed, yday, stat, macro) %>% 
  reframe(mean = sum(mean, na.rm = T),
          sd = mean * sum(sd, na.rm = T)) %>% 
  arrow::write_parquet(file.path(proc_path, "total_excr_allfarms.parquet"))

# # Quick plot
# excr %>% 
#   dplyr::filter(farm_ID == 722) %>% 
#   ggplot(aes(x = yday, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed)) +
#   geom_line(linewidth = 0.75) +
#   geom_ribbon(alpha = 0, linetype = "dotted") +
#   facet_wrap(~macro) +
#   scale_colour_manual(values = c("#666", "#cc0088", "#00aaaa")) +
#   theme_classic()

rm(excr)

## Together - difference -------------------------------------------------------------------------------------------
uneat <- file.path(proc_path, "total_uneat_allfarms.parquet") %>% arrow::read_parquet()
excr <- file.path(proc_path, "total_excr_allfarms.parquet") %>% arrow::read_parquet()
coho_biom <- file.path(proc_path, "farm_biomass_produced.qs") %>% qread()
points <- coho_biom %>% 
  distinct(farm_ID, country, geometry)

macros <- levels(uneat$macro)

uneat_diff <- uneat_diff_pc <- uneat_diff_kg.t <- uneat_diff_kg.t_pc <- list()
excr_diff <- excr_diff_pc <- excr_diff_kg.t <- excr_diff_kg.t_pc <- list()
for (m in seq_along(macros)) {
  # Mean difference in daily input rate (kg/day) between feeds across space (farms) and time (yday)
  uneat_diff[[m]] <- uneat %>% 
    dplyr::filter(macro == macros[m]) %>% 
    dplyr::select(-sd) %>% 
    mutate(mean = set_units(mean, "g") %>% set_units("kg")) %>% 
    pivot_wider(names_from = feed, values_from = mean, id_cols = c(farm_ID, yday, stat)) %>% 
    mutate(past = past-reference,
           future = future-reference) %>% 
    dplyr::select(-reference) %>% 
    pivot_longer(names_to = "feed", values_to = "change", cols = c(past, future), names_transform = list(feed = as.factor)) %>% 
    mutate(macro = macros[m])
  
  excr_diff[[m]] <- excr %>% 
    dplyr::filter(macro == macros[m]) %>% 
    dplyr::select(-sd) %>% 
    mutate(mean = set_units(mean, "g") %>% set_units("kg")) %>% 
    pivot_wider(names_from = feed, values_from = mean, id_cols = c(farm_ID, yday, stat)) %>% 
    mutate(past = past-reference,
           future = future-reference) %>% 
    dplyr::select(-reference) %>% 
    pivot_longer(names_to = "feed", values_to = "change", cols = c(past, future), names_transform = list(feed = as.factor)) %>% 
    mutate(macro = macros[m])

  # Mean percentage difference in daily input rate between feeds across space (farms) and time (yday)
  # The percentage change actually stays very consistent across time, so have condensed into a mean, min and max across time
  uneat_diff_pc[[m]] <- uneat %>% 
    dplyr::filter(macro == macros[m]) %>% 
    dplyr::select(-sd) %>% 
    pivot_wider(names_from = feed, values_from = mean, id_cols = c(farm_ID, yday, stat)) %>% 
    mutate(past = (past-reference)/reference,
           future = (future-reference)/reference) %>% 
    dplyr::select(-reference) %>% 
    pivot_longer(names_to = "feed", values_to = "perc_change", cols = c(past, future), names_transform = list(feed = as.factor)) %>% 
    group_by(farm_ID, stat, feed) %>% 
    reframe(mean_pc = mean(perc_change),
            min_pc = min(perc_change),
            max_pc = max(perc_change)) %>% 
    mutate(macro = macros[m])

  excr_diff_pc[[m]] <- excr %>% 
    dplyr::filter(macro == macros[m]) %>% 
    dplyr::select(-sd) %>% 
    pivot_wider(names_from = feed, values_from = mean, id_cols = c(farm_ID, yday, stat)) %>% 
    mutate(past = (past-reference)/reference,
           future = (future-reference)/reference) %>% 
    dplyr::select(-reference) %>% 
    pivot_longer(names_to = "feed", values_to = "perc_change", cols = c(past, future), names_transform = list(feed = as.factor)) %>% 
    group_by(farm_ID, stat, feed) %>% 
    reframe(mean_pc = mean(perc_change),
            min_pc = min(perc_change),
            max_pc = max(perc_change)) %>% 
    mutate(macro = macros[m])
  
  # Mean difference in total input per biomass produced (kg/tonnes) between feeds across space (farms) and time (yday)
  uneat_diff_kg.t[[m]] <- uneat %>% 
    dplyr::filter(macro == macros[m]) %>% 
    # Pair (cumulative) farm biomass with daily uneaten
    merge(coho_biom, by = c("farm_ID", "feed", "yday"), all = T) %>% 
    rename(input = mean.x, biomass = mean.y) %>% 
    # Add the cohorts together when they overlap
    group_by(farm_ID, feed, yday, stat, macro) %>% 
    reframe(input = sum(input),
            biomass = sum(biomass)) %>% 
    # Change uneaten units from g/day to kg/day
    mutate(input = input %>% set_units("g") %>% set_units("kg") %>% drop_units(),
           input_biom_kg.t = input/biomass) %>% 
    mutate(macro = macros[m])

  excr_diff_kg.t[[m]] <- excr %>% 
    dplyr::filter(macro == macros[m]) %>% 
    # Pair (cumulative) farm biomass with daily uneaten
    merge(coho_biom, by = c("farm_ID", "feed", "yday"), all = T) %>% 
    rename(input = mean.x, biomass = mean.y) %>% 
    # Add the cohorts together when they overlap
    group_by(farm_ID, feed, yday, stat, macro) %>% 
    reframe(input = sum(input),
            biomass = sum(biomass)) %>% 
    # Change uneaten units from g/day to kg/day
    mutate(input = input %>% set_units("g") %>% set_units("kg") %>% drop_units(),
           input_biom_kg.t = input/biomass) %>% 
    mutate(macro = macros[m])
  
  # Mean percentage difference in total input per biomass produced between feeds across space (farms) and time (yday)
  uneat_diff_kg.t_pc[[m]] <- uneat_diff_kg.t[[m]] %>% 
    # Show past and future inputs as % change from reference
    pivot_wider(names_from = feed, values_from = input_biom_kg.t, id_cols = c(farm_ID, stat, macro, yday)) %>% 
    mutate(past = (past-reference)/reference,
           future = (future-reference)/reference) %>% 
    dplyr::select(-reference) %>% 
    pivot_longer(names_to = "feed", values_to = "input_biom_pc", cols = c(past, future), names_transform = list(feed = as.factor)) %>% 
    mutate(macro = macros[m])

  excr_diff_kg.t_pc[[m]] <- excr_diff_kg.t[[m]] %>% 
    # Show past and future inputs as % change from reference
    pivot_wider(names_from = feed, values_from = input_biom_kg.t, id_cols = c(farm_ID, stat, macro, yday)) %>% 
    mutate(past = (past-reference)/reference,
           future = (future-reference)/reference) %>% 
    dplyr::select(-reference) %>% 
    pivot_longer(names_to = "feed", values_to = "input_biom_pc", cols = c(past, future), names_transform = list(feed = as.factor)) %>% 
    mutate(macro = macros[m])
}

# Put them all together
uneat_diff <- uneat_diff %>% bind_rows()
uneat_diff_pc <- uneat_diff_pc %>% bind_rows()
uneat_diff_kg.t <- uneat_diff_kg.t %>% bind_rows()
uneat_diff_kg.t_pc <- uneat_diff_kg.t_pc %>% bind_rows()
excr_diff <- excr_diff %>% bind_rows()
excr_diff_pc <- excr_diff_pc %>% bind_rows()
excr_diff_kg.t <- excr_diff_kg.t %>% bind_rows()
excr_diff_kg.t_pc <- excr_diff_kg.t_pc %>% bind_rows()

# Include geometry for plotting
uneat_diff <- uneat_diff %>% 
  mutate(macro = as.factor(macro),
         yday = as.integer(yday)) %>% 
  merge(points, by = c("farm_ID"))
uneat_diff_pc <- uneat_diff_pc %>% 
  mutate(macro = as.factor(macro)) %>% 
  merge(points, by = c("farm_ID"))
uneat_diff_kg.t <- uneat_diff_kg.t %>% 
  mutate(macro = as.factor(macro),
         yday = as.integer(yday)) %>% 
  merge(points, by = c("farm_ID"))
uneat_diff_kg.t_pc <- uneat_diff_kg.t_pc %>% 
  mutate(macro = as.factor(macro),
         yday = as.integer(yday)) %>% 
  merge(points, by = c("farm_ID"))
excr_diff <- excr_diff %>% 
  mutate(macro = as.factor(macro),
         yday = as.integer(yday)) %>% 
  merge(points, by = c("farm_ID"))
excr_diff_pc <- excr_diff_pc %>% 
  mutate(macro = as.factor(macro)) %>% 
  merge(points, by = c("farm_ID"))
excr_diff_kg.t <- excr_diff_kg.t %>% 
  mutate(macro = as.factor(macro),
         yday = as.integer(yday)) %>% 
  merge(points, by = c("farm_ID"))
excr_diff_kg.t_pc <- excr_diff_kg.t_pc %>% 
  mutate(macro = as.factor(macro),
         yday = as.integer(yday)) %>% 
  merge(points, by = c("farm_ID"))

input_diff <- rbind(uneat_diff, excr_diff)
input_diff_pc <- rbind(uneat_diff_pc, excr_diff_pc)
input_diff_kg.t <- rbind(uneat_diff_kg.t, excr_diff_kg.t)
input_diff_kg.t_pc <- rbind(uneat_diff_kg.t_pc, excr_diff_kg.t_pc)

qsave(uneat_diff, file.path(proc_path, "uneat_diff.qs"))
qsave(uneat_diff_pc, file.path(proc_path, "uneat_diff_pc.qs"))
qsave(uneat_diff_kg.t, file.path(proc_path, "uneat_diff_kg.t.qs"))
qsave(uneat_diff_kg.t_pc, file.path(proc_path, "uneat_diff_kg.t_pc.qs"))
qsave(excr_diff, file.path(proc_path, "excr_diff.qs"))
qsave(excr_diff_pc, file.path(proc_path, "excr_diff_pc.qs"))
qsave(excr_diff_kg.t, file.path(proc_path, "excr_diff_kg.t.qs"))
qsave(excr_diff_kg.t_pc, file.path(proc_path, "excr_diff_kg.t_pc.qs"))
qsave(input_diff, file.path(proc_path, "input_diff.qs"))
qsave(input_diff_pc, file.path(proc_path, "input_diff_pc.qs"))
qsave(input_diff_kg.t, file.path(proc_path, "input_diff_kg.t.qs"))
qsave(input_diff_kg.t_pc, file.path(proc_path, "input_diff_kg.t_pc.qs"))

# Quick plots
input_diff %>%
  dplyr::filter(country == "Australia") %>%
  group_by(country, macro, feed, stat, yday) %>% 
  reframe(mean = mean(change), min = min(change), max = max(change)) %>% 
  ggplot(aes(x = yday, y = mean, ymin = min, ymax = max, colour = feed, fill = feed)) +
  geom_line(linewidth = 0.75) +
  geom_ribbon(alpha = 0.2) +
  scale_colour_manual(values = c("past" = "#cc0088", "future" = "#00aaaa")) +
  scale_fill_manual(values = c("past" = "#cc0088", "future" = "#00aaaa")) +
  facet_grid(rows = vars(stat), cols = vars(macro)) +
  theme_classic() +
  theme(legend.position = "none")

uneat_diff_kg.t_pc %>%
  group_by(farm_ID) %>%
  mutate(end_day = max(yday)) %>% 
  ungroup() %>% 
  dplyr::filter(yday == end_day) %>% 
  group_by(stat, macro, feed, country) %>% 
  reframe(mean = mean(input_biom_pc),
          min = min(input_biom_pc),
          max = max(input_biom_pc)) %>% 
  ggplot(aes(x = country, y = mean, ymin = min, ymax = max, fill = feed)) +
  geom_col(position = position_dodge(), width = 0.95) +
  geom_errorbar(position = position_dodge(width = 0.95), width = 0.3) +
  scale_fill_manual(values = c("past" = "#cc0088", "future" = "#00aaaa")) +
  facet_grid(~macro) +
  theme_classic()







