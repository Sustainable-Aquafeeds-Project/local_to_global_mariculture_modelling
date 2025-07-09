# Disable all linters for this file
# nolint start

library(arrow)
library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(magrittr)
library(furrr)
library(future)
library(tictoc)
library(ggplot2)
library(fs)
library(conflicted)
library(stringr)
library(readxl)
library(units)
library(qs)
library(here)
library(targets)
conflicted::conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

# Configuration ---------------------------------------------------------------------------------------------------
# Set up parallel processing
# plan(multisession, workers = parallelly::availableCores()-1)

# Filenames
species_params_excel <- c(file = file.path(input_species_param_path, "Species.xlsx"), sheet = "Atlantic salmon")
pop_params_excel <- c(file = file.path(input_species_param_path, "Population.xlsx"))
farm_locations_parquet <- file.path(input_farm_coords_path, "farm_coords.parquet")
farm_coords_file <- file.path(output_farm_data_path, "farm_coords.qs")
farm_geometry_file <- file.path(output_farm_data_path, "farm_geometry.qs")
farm_ts_data_file <- file.path(output_farm_data_path, "farm_ts_data.qs")
species_params_file <- file.path(output_species_data_path, "species_params.qs")
sens_params_file <- file.path(output_species_data_path, "sens_params.qs")
pop_params_file <- file.path(output_species_data_path, "pop_params.qs")
feed_params_file <- file.path(output_species_data_path, "feed_params.qs")
farm_harvest_file <- file.path(output_farm_data_path, "farm_harvest_size.qs")

# Farm coordinates ------------------------------------------------------------------------------------------------
farm_coords <- if (!file.exists(farm_coords_file)) {
  times_N <- c("t_start" = 121, "t_end" = 121+547, "dt" = 1)
  times_S <- c("t_start" = 274, "t_end" = 274+547, "dt" = 1)
  
  fc <- farm_locations_parquet %>% 
    read_parquet() %>% 
    mutate(t_start = case_when(lat > 0 ~ times_N['t_start'], TRUE ~ times_S['t_start']), 
           t_end = case_when(lat > 0 ~ times_N['t_end'], TRUE ~ times_S['t_end']),
           t_start = unname(t_start),
           t_end = unname(t_end))
  
  qsave(fc, farm_coords_file)
  fc
} else {
  qread(farm_coords_file)
}

# Also save geometry for later
if (!file.exists(farm_geometry_file)) {
  file.path(input_farm_coords_path, "atlantic_salmon_locations_w_temps.qs") %>% 
    qread() %>% 
    dplyr::filter(day == "day_1") %>% 
    dplyr::select(farm_id, geometry, country) %>% 
    qsave(farm_geometry_file)
}

# Farm data -------------------------------------------------------------------------------------------------------
farm_ts_data <- if (!file.exists(farm_ts_data_file)) {
  farms_to_omit <- qread(sprintf(file.path(input_farm_coords_path, "%s_farms_to_omit.qs"), this_species))
  farm_SST_data <- read_parquet(file.path(input_farm_sst_path, "farm_SST_extracted.parquet"))
  farm_IDs <- farm_SST_data %>%
    filter(!farm_id %in% farms_to_omit) %>%
    distinct(farm_id) %>%
    pull(farm_id)
  df <- farm_SST_data %>%
    rename(farm_ID = farm_id) %>% 
    select(c(farm_ID, day, temp_c)) %>%
    mutate(day = str_split_i(day, "day_", 2) %>% as.integer())
  qsave(df, farm_ts_data_file)
  df
} else {
  qread(farm_ts_data_file)
}

# Species parameters ----------------------------------------------------------------------------------------------
species_params <- if (!file.exists(species_params_file)) {
  df <- readxl::read_excel(path = species_params_excel["file"], sheet = species_params_excel["sheet"])
  vc <- df$Value
  names(vc) <- df$Quantity
  vc[!is.na(vc)]
  qsave(vc, species_params_file)
  vc
} else {
  qread(species_params_file)
}

# Population parameters -------------------------------------------------------------------------------------------
pop_params <- if (!file.exists(pop_params_file)) {
  df <- readxl::read_excel(path = pop_params_excel["file"])
  vc <- df$Value
  names(vc) <- df$Quantity
  vc[!is.na(vc)]
  qsave(vc, pop_params_file)
  vc
} else {
  qread(pop_params_file)
}

# Feed parameters -------------------------------------------------------------------------------------------------
# Can't overwrite from here - see formulating_feeds.R
feed_params <- feed_params_file %>% qread()
reference_feed <- feed_params[["reference"]]

# Farm harvest size -----------------------------------------------------------------------------------------------
farm_harvest_size <- if (!file.exists(farm_harvest_file)) {
  farm_IDs <- farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID)
  hs <- purrr::map(farm_IDs, function(farm_id) {
    farm_times <- c(
      t_start = farm_coords$t_start[farm_coords$farm_ID == farm_id],
      t_end = farm_coords$t_end[farm_coords$farm_ID == farm_id],
      dt = 1
    )
    farm_temp <- farm_ts_data %>%
      filter(farm_ID == farm_id) %>%
      filter(day >= farm_times['t_start'], day <= farm_times['t_end']) %>%
      pull(temp_c)
    df <- fish_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = farm_temp,
      feed_params = feed_params[["reference"]],
      times = farm_times,
      init_weight = pop_params["meanW"],
      ingmax = pop_params["meanImax"]
    )
    df %>%
      as.data.frame() %>%
      slice_tail(n = 1) %>%
      select(weight) %>%
      mutate(farm_ID = farm_id,
             weight = units::set_units(weight, "g"))
  }) %>% bind_rows()
  qsave(hs, farm_harvest_file)
  hs
} else {
  qread(farm_harvest_file)
}

# Sensitivities ---------------------------------------------------------------------------------------------------
sens_all_params <- if (!file.exists(sens_params_file)) {
  vc <- c(species_params, pop_params)
  omit_names <- c('betaprot', 'betalip', 'betacarb', 'fcr', 'CS', 'nruns')
  nms <- names(vc)[!names(vc) %in% omit_names]
  vc <- vc[!names(vc) %in% omit_names]
  names(vc) <- nms
  qsave(vc, sens_params_file)
  vc
} else {
  qread(sens_params_file)
}
sens_params_names <- names(sens_all_params)
factors <- c(0.9, 1, 1.1)

Sys.setenv(TAR_PROJECT = "project_sensitivities")
rm(list = grep("tar_", ls(), value = TRUE), envir = .GlobalEnv)
targets::tar_validate()
targets::tar_outdated(callr_function = NULL)
targets::tar_make(
  # names = "tar_sens_results_spec",
  callr_function = NULL
)

overwrite <- T
mani <- tar_manifest()
sens_results <- rbind(
  tar_read(tar_sens_results_spec),
  tar_read(tar_sens_results_pop)
  )
sens_measures <- levels(sens_results$measure)
sens_results_files <- file.path(output_sens_data_path, paste0("sens_results_", sens_measures, ".qs"))
sens_results_figfiles <- file.path(output_sens_data_path, paste0("sens_plot_", sens_measures, ".qs"))
for (sm in seq_along(sens_measures)) {
  sens_results %>% 
    filter(measure == sens_measures[sm]) %>% 
    qsave(sens_results_files[sm])
  
  p <- sens_results %>% 
    filter(measure == sens_measures[sm]) %>% 
    ggplot(aes(x = adj_param, y = mean_sens, ymin = mean_sens-sd_sens, ymax = mean_sens+sd_sens)) +
    geom_col(fill = "salmon", alpha = 0.65, colour = "black") +
    geom_errorbar(width = 0.5) +
    coord_flip() +
    theme_classic()
  qsave(p, sens_results_figfiles[sm])
}

# Finish up -------------------------------------------------------------------------------------------------------
# plan(sequential)
message("All done!")

# nolint end
