# Disable all linters for this file
# nolint start

suppressMessages(suppressWarnings(suppressPackageStartupMessages({
  library(arrow)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(terra)
  library(magrittr)
  library(furrr)
  library(future)
  library(tictoc)
  library(fs)
  library(conflicted)
  library(stringr)
  library(readxl)
  library(units)
  library(qs)
  library(here)
  library(progressr)
  conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)
})))

# Configuration ---------------------------------------------------------------------------------------------------
# Set up parallel processing
# plan(multisession, workers = parallelly::availableCores()-2)
plan(multisession, workers = 8)

# Basic configuration
overwrite <- FALSE
this_species <- "atlantic_salmon"

# Directory structure
output_path <- here() %>% file.path("outputs")
gendata_path <- here() %>% file.path("data", "_general_data")
species_path <- here() %>% file.path("data", this_species)

# Input paths
input_farm_coords_path <- file.path(gendata_path, "farm_locations")
input_farm_sst_path <- file.path(gendata_path, "SST")
input_species_param_path <- file.path(species_path, "params")
input_feed_profile_path <- file.path(gendata_path, "diets")

# Output paths
output_farm_data_path <- file.path(output_path, "farm_data")
output_species_data_path <- file.path(output_path, "species_data")
output_growth_data_path <- file.path(output_path, "farm_growth_data")
output_cohorts_data_path <- file.path(output_path, "cohort_growth_data")
output_model_farm_path <- file.path(output_path, "all_outputs_farm")
output_model_cohort_path <- file.path(output_path, "all_outputs_cohort")
total_uneaten_path <- file.path(output_path, "total_uneaten_cohort")
total_excreted_path <- file.path(output_path, "total_excreted_cohort")

# Create output directories
dir_create(c(output_farm_data_path, 
             output_species_data_path, 
             output_growth_data_path,
             output_cohorts_data_path,
             output_model_farm_path,
             output_model_cohort_path,
             total_uneaten_path,
             total_excreted_path))

# Print configuration summary
message(sprintf("Configuration:\n- Species: %s\n- Parallel workers: %d\n- Overwrite existing: %s\n- Output path: %s",
               this_species,
               works,
               ifelse(overwrite, "yes", "no"),
               output_path))

# Helper functions ------------------------------------------------------------------------------------------------
file_exists_skip <- function(filepath, section_name, fn) {
  if (file.exists(filepath) && !overwrite) {
    message(sprintf("Skipping %s - file exists", section_name))
    return(qs::qread(filepath))
  } else {
    tic()
    result <- fn()
    qs::qsave(result, filepath)
    message(sprintf("Saved %s - %s", section_name, toc(quiet = T)$callback_msg))
    return(result)
  }
}

source("00_model_functions.R")

# Load farm harvest size data -------------------------------------------------------------------------------------
# This is from the individual running
farm_harvest_size <- qs::qread(file.path(output_farm_data_path, "farm_harvest_size.qs")) %>%
  select(c(farm_ID, weight)) %>%
  mutate(weight = units::set_units(weight, "g"))

# Farm coordinates ------------------------------------------------------------------------------------------------
farm_coords_file <- file.path(output_farm_data_path, "farm_coords.qs")
farm_coords <- file_exists_skip(farm_coords_file, "farm coordinates", function() {
    times_N <- c("t_start" = 121, "t_end" = 121+547, "dt" = 1)
    times_S <- c("t_start" = 274, "t_end" = 274+547, "dt" = 1)
    
    read_parquet(file.path(input_farm_coords_path, "farm_coords.parquet")) %>% 
      mutate(t_start = case_when(lat > 0 ~ times_N['t_start'], TRUE ~ times_S['t_start']), 
             t_end = case_when(lat > 0 ~ times_N['t_end'], TRUE ~ times_S['t_end']),
             t_start = unname(t_start),
             t_end = unname(t_end))
  })

# Also save geometry for later
farm_geometry_file <- file.path(output_farm_data_path, "farm_geometry.qs")
farm_coords <- file_exists_skip(farm_geometry_file, "farm geometry", function() {
  file.path(input_farm_coords_path, "atlantic_salmon_locations_w_temps.qs") %>% 
    qread() %>% 
    dplyr::filter(day == "day_1") %>% 
    dplyr::select(farm_id, geometry, country)
  })

# Farm static data ------------------------------------------------------------------------------------------------
farm_static_data_file <- file.path(output_farm_data_path, "farm_static_data.qs")
farm_static_data <- file_exists_skip(farm_static_data_file, "farm static data", function() {
    farms_to_omit <- qs::qread(sprintf(file.path(input_farm_coords_path, "%s_farms_to_omit.qs"), this_species))
    farm_SST_data <- read_parquet(file.path(input_farm_sst_path, "farm_SST_extracted.parquet"))
    farm_IDs <- farm_SST_data %>%
      filter(!farm_id %in% farms_to_omit) %>%
      distinct(farm_id) %>%
      pull(farm_id)
    
    farm_SST_data %>%
      rename(farm_ID = farm_id) %>% 
      distinct(iso3c, country, tonnes_per_farm, daily_mort_rate, farm_ID) %>%
      mutate(tonnes_per_farm = tonnes_per_farm %>% units::set_units("t") %>% units::set_units("g"))
  })

# Farm timeseries data --------------------------------------------------------------------------------------------
farm_ts_data_file <- file.path(output_farm_data_path, "farm_ts_data.qs")
farm_ts_data <- file_exists_skip(farm_ts_data_file, "farm ts data", function() {
    farms_to_omit <- qs::qread(sprintf(file.path(input_farm_coords_path, "%s_farms_to_omit.qs"), this_species))
    
    farm_SST_data <- read_parquet(file.path(input_farm_sst_path, "farm_SST_extracted.parquet"))
    farm_IDs <- farm_SST_data %>%
      filter(!farm_id %in% farms_to_omit) %>%
      distinct(farm_id) %>%
      pull(farm_id)
    
    farm_SST_data %>%
      rename(farm_ID = farm_id) %>% 
      select(c(farm_ID, day, temp_c)) %>%
      mutate(day = str_split_i(day, "day_", 2) %>% as.numeric())
  })

# Species parameters ----------------------------------------------------------------------------------------------
species_params_file <- file.path(output_species_data_path, "species_params.qs")
species_params <- file_exists_skip(species_params_file, "species parameters", function() {
    df <- readxl::read_excel(file.path(input_species_param_path, "Species.xlsx"), sheet = "Atlantic salmon")
    vc <- df$Value
    names(vc) <- df$Quantity
    vc[!is.na(vc)]
  })

# Population parameters -------------------------------------------------------------------------------------------
pop_params_file <- file.path(output_species_data_path, "pop_params.qs")
if (!exists("pop_params")) {
  pop_params <- if (file.exists(pop_params_file)) {
    qs::qread(pop_params_file)
  } else {
    file_exists_skip(pop_params_file, "population parameters", function() {
      df <- readxl::read_excel(file.path(input_species_param_path, "Population.xlsx"))
      vc <- df$Value
      names(vc) <- df$Quantity
      vc[!is.na(vc)]
    })
  }
}

# Feed parameters -------------------------------------------------------------------------------------------------
feed_params_file <- file.path(output_species_data_path, "feed_params.qs")
if (!exists("feed_params")) {
  feed_params <- if (file.exists(feed_params_file)) {
    qs::qread(feed_params_file)
  } else {
    file_exists_skip(feed_params_file, "feed parameters", function() {
      feed_types <- c("reference", "past", "future")
      feed_profile_file <- file.path(input_feed_profile_path, "feed_profiles_and_composition_Atlantic_salmon.xlsx")
      
      future_map(feed_types, function(ft) {
        df <- readxl::read_excel(feed_profile_file, sheet = ft)
        list(
          Proteins = df %>% select(ingredient, proportion, contains("protein")) %>%
            select(-contains("feed")) %>%
            rename(macro = ing_protein, digest = ing_protein_digestibility),
          Carbohydrates = df %>% select(ingredient, proportion, contains("carb")) %>%
            select(-contains("feed")) %>%
            rename(macro = ing_carb, digest = ing_carb_digestibility),
          Lipids = df %>% select(ingredient, proportion, contains("lipid")) %>%
            select(-contains("feed")) %>%
            rename(macro = ing_lipid, digest = ing_lipid_digestibility)
        )
      }) %>% setNames(feed_types)
    })
  }
}

# Main farm growth ------------------------------------------------------------------------------------------------
farm_IDs <- farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID)
farm_growth_files <- file.path(output_growth_data_path, sprintf("main_farm_growth_%s.qs", fixnum(farm_IDs)))

message("Processing main farm growth...")
tic("Main farm growth")
with_progress({
  p <- progressor(along = farm_growth_files)
  future_map2(farm_IDs, farm_growth_files, function(farm_id, output_file) {
    if (!file.exists(output_file) || overwrite) {
    result <- future_map(feed_params, function(feed_type) {
      farm_times <- c(
        t_start = farm_coords$t_start[farm_coords$farm_ID == farm_id],
        t_end = farm_coords$t_end[farm_coords$farm_ID == farm_id],
        dt = 1
      )
      
      farm_temp <- farm_ts_data %>%
        filter(farm_ID == farm_id) %>%
        filter(day >= farm_times['t_start'], day <= farm_times['t_end']) %>%
        pull(temp_c)
      
      N_pop <- generate_pop(
        harvest_n = 1.117811 * units::drop_units(
          farm_static_data$tonnes_per_farm[farm_static_data$farm_ID == farm_id] /
            farm_harvest_size$weight[farm_harvest_size$farm_ID == farm_id]
        ),
        mort = pop_params['mortmyt'],
        times = farm_times
      )
      
      farm_growth(
        pop_params = pop_params,
        species_params = species_params,
        water_temp = farm_temp,
        feed_params = feed_type,
        times = farm_times,
        N_pop = N_pop,
        nruns = pop_params['nruns']
      )
    }, .options = furrr_options(seed = TRUE)) %>% 
      setNames(names(feed_params))
    
    qs::qsave(result, output_file)
    }
    p(sprintf("Saved %s", basename(output_file)))
  }, .options = furrr_options(seed = TRUE))
})
toc()

# Cohorts data ----------------------------------------------------------------------------------------------------
cohorts_files <- file.path(output_cohorts_data_path, sprintf("cohorts_data_%s.qs", fixnum(farm_IDs)))
stat_names <- qs::qread(farm_growth_files[1])[[1]] %>% names()

tic("Cohort growth")
with_progress({
  p <- progressor(along = farm_growth_files)
  future_map2(farm_growth_files, cohorts_files, function(farm_file, cohorts_file) {
  if (all(!file.exists(cohorts_file), file.exists(farm_file), !overwrite)) {
    file_exists_skip(cohorts_file, "cohorts data", function() {
      farm_data <- qs::qread(farm_file)
      farm_id <- farm_file %>% str_split_i("[.]", 1) %>% str_split_i("main_farm_growth_", 2) %>% as.integer()
      
      farm_times <- c(
        t_start = farm_coords$t_start[farm_coords$farm_ID == farm_id],
        t_end = farm_coords$t_end[farm_coords$farm_ID == farm_id],
        dt = 1
      )
      
      ls <- future_map(stat_names, function(stat) {
        df0 <- cbind(matrix(farm_times["t_start"]:farm_times["t_end"], 548, 1), matrix(1, 548, 1))
        df0 <- cbind(df0, matrix(1:nrow(df0), 548, 1))
        df1 <- df2 <- df3 <- cbind(df0, farm_data[[feed]][[stat]])
        colnames(df1) <- colnames(df2) <- colnames(df3) <- c("yday", "cohort", "prod_day", "mean", "sd")
        df2[,'cohort'] <- 2
        df2[,'yday'] <- df2[,'yday']+365
        df3[,'cohort'] <- 3
        df3[,'yday'] <- df3[,'yday']+365+365
        
        rbind(df1, df2, df3)
      }, .options = furrr_options(seed = TRUE)) 
      

        setNames(names(stat_names))
    })
  }
}, .progress = TRUE, .options = furrr_options(seed = TRUE))
})
toc()

# Feed differences ------------------------------------------------------------------------------------------------
message("Processing feed differences...")
cohort_data <- cohorts_files[1] %>% qs::qread()

tic("Model outputs processing")
with_progress({
  p <- progressor(along = farm_growth_files)
  future_map(seq_along(farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID)), function(fid) {
  farm_id <- farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID) %>% .[[fid]]
  
  # Process farm-level statistics
  fnames_farms <- file.path(output_model_farm_path, 
                           paste0(paste("farmID", fixnum(farm_id), farm_stats, sep = "_"), ".parquet"))
  
  if (overwrite || any(!file.exists(fnames_farms))) {
    walk(seq_along(farm_stats), function(st) {
      if (overwrite || !file.exists(fnames_farms[st])) {
        map_dfr(names(feed_params), function(feed_type) {
          main_farm_growth[[feed_type]][[fid]][[farm_stats[st]]] %>%
            as.data.frame() %>%
            rename(mean = V1, sd = V2) %>%
            mutate(feed = feed_type, days = 1:548)
        }) %>%
          relocate(days, .before = mean) %>%
          mutate(farm_ID = farm_id,
                feed = factor(feed, levels = names(feed_params))) %>%
          write_parquet(fnames_farms[st])
      }
    })
  }
  
  # Process cohort statistics
  fnames_cohorts <- file.path(output_model_cohort_path, 
                             paste0(paste("farmID", fixnum(farm_id), "cohorts", cohort_stats, sep = "_"), ".parquet"))
  
  if (overwrite || any(!file.exists(fnames_cohorts))) {
    walk(seq_along(cohort_stats), function(st) {
      if (overwrite || !file.exists(fnames_cohorts[st])) {
        map_dfr(names(feed_params), function(feed_type) {
          filter(cohorts_data[[feed_type]][[fid]], stat == cohort_stats[st])
        }) %>%
          mutate(feed = factor(feed, levels = names(feed_params))) %>%
          write_parquet(fnames_cohorts[st])
      }
    })
  }
}, .progress = TRUE)
})
toc()

# Total uneaten data ----------------------------------------------------------------------------------------------
message("Processing total uneaten data...")
tic("Total uneaten processing")
with_progress({
  p <- progressor(along = farm_growth_files)
  future_map(seq_along(farm_data$farm_IDs), function(fid) {
  farm_id <- farm_data$farm_IDs[fid]
  fname <- file.path(total_uneaten_path, paste0("farmID_", fixnum(farm_id), "_total_uneat.parquet"))
  
  if (overwrite || !file.exists(fname)) {
    map_dfr(names(feed_params), function(feed_type) {
      cohorts_data[[feed_type]][[fid]] %>%
        filter(str_detect(stat, "_uneat$")) %>%
        mutate(macro = str_split_i(stat, "_", 1),
               stat = str_split_i(stat, "_", 2)) %>%
        mutate(macro = as.factor(macro),
               stat = factor(stat, levels = c("uneat", "excr"))) %>%
        write_parquet(fname)
    })
  }
}, .progress = TRUE)
})
toc()

# Total excreted data ---------------------------------------------------------------------------------------------
message("Processing total excreted data...")
tic("Total excreted processing")
with_progress({
  p <- progressor(along = farm_growth_files)
  future_map(seq_along(farm_data$farm_IDs), function(fid) {
  farm_id <- farm_data$farm_IDs[fid]
  fname <- file.path(total_excreted_path, paste0("farmID_", fixnum(farm_id), "_total_excr.parquet"))
  
  if (overwrite || !file.exists(fname)) {
    map_dfr(names(feed_params), function(feed_type) {
      cohorts_data[[feed_type]][[fid]] %>%
        filter(str_detect(stat, "_excr$")) %>%
        mutate(macro = str_split_i(stat, "_", 1),
               stat = str_split_i(stat, "_", 2)) %>%
        mutate(macro = as.factor(macro),
               stat = factor(stat, levels = c("uneat", "excr"))) %>%
        write_parquet(fname)
    })
  }
}, .progress = TRUE)
})
toc()

message("Completed successfully!")

# nolint end
