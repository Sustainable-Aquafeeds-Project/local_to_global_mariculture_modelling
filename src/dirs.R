# This sets up all the directories for the project so that all scripts use the same directory names
suppressPackageStartupMessages(suppressWarnings({
library(here)
library(stringr)
}))

# From run_individual
this_species <- "atlantic_salmon"

# Directory structure
output_path <- here() %>% file.path("outputs")
gendata_path <- here() %>% file.path("data", "_general_data")
rawdata_path <- here() %>% file.path("data", "raw")
bigdata_path <- "D:/local_to_global_mariculture_modelling_data"
rdsi_dir <- "/mnt/rdsi"
species_path <- here() %>% file.path("data", this_species)
species_fig_path <- file.path(species_path, "figures")

# Input paths
input_farm_coords_path <- file.path(gendata_path, "farm_locations")
input_farm_sst_path <- file.path(gendata_path, "SST")
input_species_param_path <- file.path(species_path, "params")
input_feed_profile_path <- file.path(gendata_path, "diets")
input_spec_layers_path <- file.path(gendata_path, "species_layers")
input_aquamaps_path <- file.path(input_spec_layers_path, "from_aquamaps")
input_Ndata_path <- file.path(gendata_path, "background_nitrogen")

# Output paths
output_farm_data_path <- file.path(output_path, "farm_data")
output_species_data_path <- file.path(output_path, "species_data")
output_sens_data_path <- file.path(output_path, "sensitivity_data")
per_biomass_farm_path <- file.path(output_path, "per_biomass_outputs_farm")
per_biomass_cohort_path <- file.path(output_path, "per_biomass_outputs_cohort")
output_model_farm_path <- file.path(output_path, "all_outputs_farm")
output_model_cohort_path <- file.path(output_path, "all_outputs_cohort")
data_analysis_path <- file.path(output_path, "data_analysis")
comparisons_path <- file.path(data_analysis_path, "feed_comparisons")
impacts_path <- file.path(output_path, "nutrient_impacts")

# Create output directories
dir.create(output_species_data_path, showWarnings = F)
dir.create(output_sens_data_path, showWarnings = F)
dir.create(output_farm_data_path, showWarnings = F)
dir.create(per_biomass_farm_path, showWarnings = F)
dir.create(per_biomass_cohort_path, showWarnings = F)
dir.create(output_model_farm_path, showWarnings = F)
dir.create(output_model_cohort_path, showWarnings = F)
dir.create(data_analysis_path, showWarnings = F)
dir.create(comparisons_path, showWarnings = F)
dir.create(impacts_path, showWarnings = F)
