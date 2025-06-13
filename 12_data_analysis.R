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

# This script aims to take the model outputs from the targets pipeline and answer specific questions.
# Paths & globals -------------------------------------------------------------------------------------------------
source("00_model_functions.R")
source("00_dirs.R")

cohort_comp_fnms <- output_cohorts_data_path %>% 
  list.files(full.names = T) %>% 
  str_subset("farmrun_comparisons")

aggregate_comparison_files <- output_cohorts_data_path %>% 
  list.files(full.names = T) %>% 
  str_subset("allfarms_comparisons")

# * Main point is the difference between feeds. 
# * This should be expressed in t/t salmon, so farm biomass is also important. 
# * Also interesting to see if difference between feeds varies geographically (correlates with mean/median/max temperature?)
remove_unit("g_fish")
remove_unit("kg_fish")
remove_unit("t_fish")
install_unit(symbol = "g_fish")
install_unit(symbol = "kg_fish", def = "1000 g_fish")
install_unit(symbol = "t_fish", def = "1000 kg_fish")

# Different feeds -------------------------------------------------------------------------------------------------
## Question 1 -----------------------------------------------------------------------------------------------------
# How does the change in feed affect daily and overall excretion and uneaten feed waste (P, L, C, total)?

coho_biom <- aggregate_comparison_files %>% 
  str_subset("biomass_stat") %>% 
  qs::qread() %>% 
  mutate(t = as.integer(t))

all_inputs <- rbind(
  aggregate_comparison_files %>% 
    str_subset("total_excr_stat") %>% qread() %>% mutate(measure = "excr", category = "total"),
  aggregate_comparison_files %>% 
    str_subset("P_excr_stat") %>% qread() %>% mutate(measure = "excr", category = "P"),
  aggregate_comparison_files %>% 
    str_subset("L_excr_stat") %>% qread() %>% mutate(measure = "excr", category = "L"),
  aggregate_comparison_files %>% 
    str_subset("C_excr_stat") %>% qread() %>% mutate(measure = "excr", category = "C"),
  aggregate_comparison_files %>% 
    str_subset("total_uneat_stat") %>% qread() %>% mutate(measure = "uneat", category = "total"),
  aggregate_comparison_files %>% 
    str_subset("P_uneat_stat") %>% qread() %>% mutate(measure = "uneat", category = "P"),
  aggregate_comparison_files %>% 
    str_subset("L_uneat_stat") %>% qread() %>% mutate(measure = "uneat", category = "L"),
  aggregate_comparison_files %>% 
    str_subset("C_uneat_stat") %>% qread() %>% mutate(measure = "uneat", category = "C")
) %>% 
  mutate(measure = as.factor(measure),
         category = as.factor(category), 
         t = as.integer(t))

all_inputs <- all_inputs %>% 
  merge(coho_biom, by = c("farm_ID", "feed", "t"), all = T) %>% 
  rename(mean = mean.x, sd = sd.x, biomass_mean = mean.y, biomass_sd = sd.y) %>% 
  mutate(mean = set_units(mean, "g d-1"),
         sd = set_units(sd, "g d-1"),
         biomass_mean = set_units(biomass_mean, "g_fish"),
         biomass_sd = set_units(biomass_sd, "g_fish"))

qs::qsave(all_inputs, file.path(data_analysis_path, "all_waste_inputs.qs"))

### Total mass lost -----------------------------------------------------------------------------------------------
total_inputs <- all_inputs %>% 
  filter(category == "total") %>% 
  group_by(farm_ID, feed, measure) %>% 
  reframe(biomass_total = max(biomass_mean),
          total = sum(mean) %>% drop_units() %>% set_units("g")) 

qs::qsave(total_inputs, file.path(data_analysis_path, "sumtotal_waste_inputs.qs"))

# Quick glance
# total_inputs %>% 
#   group_by(feed, measure) %>% 
#   reframe(protein_lost = mean(P_total/biomass_total)) %>% 
#   mutate(protein_lost = protein_lost %>% set_units("g kg_fish-1")) %>% 
#   pivot_wider(names_from = measure, values_from = protein_lost) %>% 
#   mutate(from_faeces = excr/(excr + uneat))

### Total protein lost --------------------------------------------------------------------------------------------
total_P_inputs <- all_inputs %>% 
  filter(category == "P") %>% 
  group_by(farm_ID, feed, measure) %>% 
  reframe(biomass_total = max(biomass_mean),
          P_total = sum(mean) %>% drop_units() %>% set_units("g")) 

qs::qsave(total_P_inputs, file.path(data_analysis_path, "sumtotal_waste_P_inputs.qs"))

# Quick glance
# total_inputs %>% 
#   group_by(feed, measure) %>% 
#   reframe(protein_lost = mean(P_total/biomass_total)) %>% 
#   mutate(protein_lost = protein_lost %>% set_units("g kg_fish-1")) %>% 
#   pivot_wider(names_from = measure, values_from = protein_lost) %>% 
#   mutate(from_faeces = excr/(excr + uneat))

## Question 2 -----------------------------------------------------------------------------------------------------
# How does the change in feed affect daily and overall N inputs (total, % from faeces)?
N_inputs <- all_inputs %>% 
  filter(category == "P") %>% 
  mutate(mean = mean * 0.16,
         sd = sd * 0.16)

qs::qsave(N_inputs, file.path(data_analysis_path, "all_N_inputs.qs"))

total_N_inputs <- N_inputs %>% 
  group_by(farm_ID, feed, measure) %>% 
  reframe(biomass_total = max(biomass_mean),
          N_total = sum(mean) %>% drop_units() %>% set_units("g")) 

qs::qsave(total_N_inputs, file.path(data_analysis_path, "sumtotal_N_inputs.qs"))

# total_N_inputs %>%
#   group_by(feed, measure) %>%
#   reframe(N_lost = mean(N_total/biomass_total)) %>%
#   mutate(N_lost = N_lost %>% set_units("g kg_fish-1")) %>%
#   pivot_wider(names_from = measure, values_from = N_lost) %>%
#   mutate(from_faeces = excr/(excr + uneat))


## Question 3 -----------------------------------------------------------------------------------------------------
# Does total_excr (g/gfish/day) change over time for any farm, and does the feed make a difference in that?

total_excr <- purrr::map_dfr(cohort_comp_fnms, function(fnm) {
  excr <- aggregate_comparison_files %>% str_subset("total_excr_stat") %>% qread()
  biom <- aggregate_comparison_files %>% str_subset("biomass_stat") %>% qread()
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


# Biomass produced ------------------------------------------------------------------------------------------------
# How accurately is the model producing the correct biomass production levels?
coho_biom <- aggregate_comparison_files %>% 
  str_subset("biomass_stat") %>% 
  qs::qread() %>% 
  filter(feed == "reference") %>% 
  group_by(farm_ID) %>% 
  slice_max(t) %>% 
  mutate(mean = mean %>% set_units("g") %>% set_units("t") %>% drop_units()) %>% 
  merge(farm_info, by.x = "farm_ID", by.y = "farm_id")

farm_info <- file.path(input_farm_coords_path, "atlantic_salmon_locations_w_temps.qs") %>% 
  qread() %>% 
  filter(day == "day_1") %>% 
  mutate(country = as.factor(country),
         day = str_split_i(day, "_", 2) %>% as.integer()) %>% 
  dplyr::select(-c(model_name, F_CODE, data_year, data_source, details, data_type_2, data_type, species_group, 
                   harvest_n, stocking_n, harvest_size_t, daily_mort_rate, day))
hist(farm_info$tonnes_per_farm)



qsave(coho_biom, file.path(data_analysis_path, "biomass_produced_comparison.qs"))

ggplot(coho_biom, aes(x = tonnes_per_farm, y = mean, colour = country)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_classic() +
  # scale_x_continuous(limits = c(0,1700)) +
  # scale_y_continuous(limits = c(0,1700)) +
  labs(y = "Modelled biomass produced (t)", x = "Observed biomass produced (t)")

ggsave(file.path(data_analysis_path, "biomass_produced_comparison_plot.png"))




