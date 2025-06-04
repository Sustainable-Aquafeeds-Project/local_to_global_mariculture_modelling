# nolint start

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
  library(units)
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

input_coho_path <- here() %>% file.path("outputs", "cohort_growth_data")
input_farm_path <- here() %>% file.path("outputs", "farm_growth_data")
cohort_comp_fnms <- input_coho_path %>% 
  list.files(full.names = T) %>% 
  str_subset("farmrun_comparisons")
cohort_refe_fnms <- input_farm_path %>% 
  list.files(full.names = T) %>% 
  str_subset("farmrun_reference")
cohort_past_fnms <- input_farm_path %>% 
  list.files(full.names = T) %>% 
  str_subset("farmrun_past")
cohort_futu_fnms <- input_farm_path %>% 
  list.files(full.names = T) %>% 
  str_subset("farmrun_future")

output_coho_path <- file.path("data", "atlantic_salmon", "data_products", "model_outputs_cohort")
data_analysis_path <- file.path("outputs", "data_analysis")
dir.create(data_analysis_path)

# * Main point is the difference between feeds. 
# * This should be expressed in t/t salmon, so farm biomass is also important. 
# * Also interesting to see if difference between feeds varies geographically (correlates with mean/median/max temperature?)
remove_unit("g_fish")
remove_unit("kg_fish")
remove_unit("t_fish")
install_unit(symbol = "g_fish", name = "grams of salmon biomass")
install_unit(symbol = "kg_fish", def = "1000 g_fish")
install_unit(symbol = "t_fish", def = "1000 kg_fish")

# Biomass produced ------------------------------------------------------------------------------------------------
# How accurately is the model producing the correct biomass production levels?
farm_info <- file.path("data", "_general_data", "farm_locations", "atlantic_salmon_locations_w_temps.qs") %>% 
  qread() %>% 
  filter(day == "day_1") %>% 
  mutate(country = as.factor(country),
         day = str_split_i(day, "_", 2) %>% as.integer()) %>% 
  dplyr::select(-c(model_name, F_CODE, data_year, data_source, details, data_type_2, data_type, species_group, 
                   harvest_n, stocking_n, harvest_size_t, daily_mort_rate, day))
hist(farm_info$tonnes_per_farm)

coho_biom <- purrr::map_dfr(cohort_refe_fnms, function(fnm) {
    qs::qread(fnm)[["biomass_stat"]] %>% as.data.frame()
  }) %>% 
  group_by(farm_ID) %>% 
  slice_max(t) %>% 
  mutate(mean = mean %>% set_units("g") %>% set_units("t") %>% drop_units()) %>% 
  merge(farm_info, by.x = "farm_ID", by.y = "farm_id")

qsave(coho_biom, file.path(data_analysis_path, "biomass_produced_comparison.qs"))

ggplot(coho_biom, aes(x = tonnes_per_farm, y = mean, colour = country)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic() +
  # scale_x_continuous(limits = c(0,1700)) +
  # scale_y_continuous(limits = c(0,1700)) +
  labs(y = "Modelled biomass produced (t)", x = "Observed biomass produced (t)")

ggsave(file.path(data_analysis_path, "biomass_produced_comparison_plot.png"))

# Different feeds -------------------------------------------------------------------------------------------------
## Question 1 -----------------------------------------------------------------------------------------------------
# Does total_excr (g/gfish/day) change over time for any farm, and does the feed make a difference in that?
total_excr <- purrr::map_dfr(cohort_comp_fnms, function(fnm) {
  excr <- qs::qread(fnm)[["total_excr_stat"]]
  biom <- qs::qread(fnm)[["biomass_stat"]]
  merge(excr, biom, by = c("farm_ID", "feed", "t")) %>% 
    rename(total_excr = mean.x, biomass = mean.y) %>% 
    dplyr::select(-c(sd.x, sd.y))
}) %>% 
  mutate(total_excr = set_units(total_excr, "g d-1"),
         biomass = set_units(biomass, "g_fish"))

over_time <- total_excr %>% 
  mutate(excr_biom = total_excr/biomass) %>% 
  group_by(farm_ID, feed) %>% 
  mutate(diff = c(NA, diff(excr_biom))) %>% 
  ungroup()

minna(over_time$diff)
maxna(over_time$diff)

over_time %>% 
  ggplot(aes(x = diff, fill = feed)) +
  geom_histogram(alpha = 0.5, position = "identity", colour = "black")

# The answer is barely
rm(total_excr, over_time)

## Question 2 -----------------------------------------------------------------------------------------------------
# How does the change in feed affect:
# * Daily and overall excretion (P, L, C, total)
# * Daily and overall uneaten feed (P, L, C, total)
# * Daily and overall inputs (excretion + uneaten feed)

# All times
all_inpts <- purrr::map_dfr(cohort_comp_fnms, function(fnm) {
  excr <- rbind(
    qs::qread(fnm)[["total_excr_stat"]] %>% mutate(measure = "excr", category = "total"),
    qs::qread(fnm)[["P_excr_stat"]] %>% mutate(measure = "excr", category = "P"),
    qs::qread(fnm)[["L_excr_stat"]] %>% mutate(measure = "excr", category = "L"),
    qs::qread(fnm)[["C_excr_stat"]] %>% mutate(measure = "excr", category = "C"),
    qs::qread(fnm)[["total_uneat_stat"]] %>% mutate(measure = "uneat", category = "total"),
    qs::qread(fnm)[["P_uneat_stat"]] %>% mutate(measure = "uneat", category = "P"),
    qs::qread(fnm)[["L_uneat_stat"]] %>% mutate(measure = "uneat", category = "L"),
    qs::qread(fnm)[["C_uneat_stat"]] %>% mutate(measure = "uneat", category = "C")
  ) %>% 
    merge(qs::qread(fnm)[["biomass_stat"]], by = c("farm_ID", "feed", "t"), all = T) %>% 
    rename(mean = mean.x, sd = sd.x, biomass_mean = mean.y, biomass_sd = sd.y) %>% 
    mutate(measure = as.factor(measure), category = as.factor(category))
}) %>% 
  mutate(mean = set_units(mean, "g d-1"),
         sd = set_units(sd, "g d-1"),
         biomass_mean = set_units(biomass_mean, "g_fish"),
         biomass_sd = set_units(biomass_sd, "g_fish"),
         mean_biom = mean/biomass_mean)

qs::qsave(all_inpts, file.path(data_analysis_path, "all_inputs.qs"))

all_inpts %>% 
  group_by(feed, t, measure, category) %>% 
  reframe(mean = meanna(mean)) %>% 
  ggplot(aes(x = t, y = mean, colour = feed)) +
  geom_line(linewidth = 0.75) +
  facet_grid(rows = vars(category), cols = vars(measure)) +
  theme_classic() +
  labs(x = "Day of the year", y = "Mean")

# Stats throughout 1 cohort production period
inpts_stats <- all_inpts %>% 
  group_by(farm_ID, feed, measure, category) %>% 
  reframe(mean_value = meanna(mean),
          sd_value = sdna(mean) %>% set_units("g d-1"),
          min_value = minna(mean),
          max_value = maxna(mean),
          mean_value_biom = meanna(mean_biom),
          sd_value_biom = sdna(mean_biom) %>% set_units("g g_fish-1 d-1"),
          min_value_biom = minna(mean_biom),
          max_value_biom = maxna(mean_biom)) %>% 
  merge(farm_info, by.x = "farm_ID", by.y = "farm_id")

qs::qsave(inpts_stats, file.path(data_analysis_path, "inputs_stats.qs"))

# Total throughout 1 cohort production period
inpts_end <- all_inpts %>% 
  group_by(farm_ID, feed, measure, category) %>% 
  reframe(sum_value = meanna(mean),
          sum_value_biom = meanna(mean_biom)) %>% 
  merge(farm_info, by.x = "farm_ID", by.y = "farm_id")

qs::qsave(inpts_end, file.path(data_analysis_path, "inputs_endsums.qs"))

