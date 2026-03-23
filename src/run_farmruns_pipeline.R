# ==========================================================================================================================================================
# run_farmruns_pipeline.R
# 
# Standalone (non-targets) version of the farm runs pipeline.
# Reproduces all outputs of targets/_targets_farmruns.R, saving each result
# as an individual .qs file keyed by farm_ID and feed_name.
#
# Output sections:
#   1. farm_runs_full  – per-fish daily timeseries (all farms × all feeds)
#   2. fish_ends       – per-fish harvest summaries (all farms × all feeds)
#   3. farm_runs       – population-weighted daily timeseries (sample farms only)
#   4. farm_ends       – population-weighted harvest summaries (all farms × all feeds)
# ==========================================================================================================================================================

# ==========================================================================================================================================================
# 0. OPTIONS
# ==========================================================================================================================================================

inds_per_farm <- 2000
farm_bite_size <- 2750

# ==========================================================================================================================================================
# 0. SETUP: packages, functions, paths
# ==========================================================================================================================================================
suppressPackageStartupMessages(suppressWarnings({
  library(qs)
  library(tidyverse)
  library(matrixStats)
  library(sf)
  library(future)
  library(furrr)
  library(tictoc)
  library(here)
}))

source(here("src/dirs.R"))
source(here("src/functions.R"))
source(here("src/model_functions_new.R"))

# Output root and sub-directories
out_root <- file.path(output_path, "model_run_outputs")

dir_farm_runs_full <- file.path(out_root, "farm_runs_full")
dir_fish_ends      <- file.path(out_root, "fish_ends")
dir_farm_runs      <- file.path(out_root, "farm_runs")
dir_farm_ends      <- file.path(out_root, "farm_ends")

bigdir_farm_runs_full <- file.path(bigdata_path, "outputs", "model_run_outputs", "farm_runs_full")

for (d in c(out_root, dir_farm_runs_full, dir_fish_ends, dir_farm_runs, dir_farm_ends)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Sometimes files get orphaned - delete all that are the wrong size
# message("\nCleaning up rogue files...")
# files <- file.info(list.files(out_root, recursive = T, full.names = TRUE, pattern = ".qs")) %>% 
#   tibble::rownames_to_column("filepath") %>% 
#   mutate(
#     filename = basename(filepath),
#     size_kb   = size / 1024,
#     min_size = case_when(
#       str_detect(filepath, dir_farm_runs_full) ~ 67000,
#       str_detect(filepath, dir_fish_ends) ~ 230,
#       str_detect(filepath, dir_farm_runs) ~ 130,
#       str_detect(filepath, dir_farm_ends) ~ 0,
#     ),
#     max_size = case_when(
#       str_detect(filepath, dir_farm_runs_full) ~ 80000,
#       str_detect(filepath, dir_fish_ends) ~ 255,
#       str_detect(filepath, dir_farm_runs) ~ 160,
#       str_detect(filepath, dir_farm_ends) ~ 10,
#     )
#   ) %>% 
#   select(filepath, filename, size_kb, min_size, max_size)

# files_to_delete <- files %>% 
#   filter(size_kb < min_size | size_kb > max_size)
# files_to_delete %>% 
#   select(filename, size_kb) %>% 
#   print()

# file.remove(files_to_delete$filepath)
# rm(files_to_delete, files)
# gc()

# ==========================================================================================================================================================
# 1. LOAD DATA
# ==========================================================================================================================================================

message("\nLoading data...")

all_params <- c(
  qread(file.path(output_species_data_path, "species_params.qs")),
  qread(file.path(output_species_data_path, "pop_params.qs"))
)

feed_params <- c(
  qread(file.path(output_species_data_path, "feed_params_MD.qs")),
  qread(file.path(output_species_data_path, "feed_params_PD.qs")),
  qread(file.path(output_species_data_path, "feed_params_AD.qs"))
)
feed_names <- names(feed_params)

# Create seperate dirs
for (d in 1:length(feed_names)) {
  dir.create(file.path(dir_farm_runs_full, feed_names[d]), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir_fish_ends, feed_names[d]), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir_farm_runs, feed_names[d]), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir_farm_ends, feed_names[d]), showWarnings = FALSE, recursive = TRUE)
}

stats_measures <- c(
  "food_prov", "total_excr", "total_uneat", "nitrogen_excr", "nitrogen_uneat", "total_nitrogen", "carbon_excr", "carbon_uneat", "total_carbon"
)

farm_static_data <- qread(file.path(input_farm_coords_path, "farm_static_data_w_stocking.qs")) %>%
  st_drop_geometry() %>%
  rename(farm_ID = farm_id)

farm_temporal_data_file <- file.path(out_root, "farm_temporal_data_withpop.qs")
if (file.exists(farm_temporal_data_file)) {
  farm_temporal_data <- qread(farm_temporal_data_file)
} else {
  farm_temporal_data_raw <- qread(file.path(input_farm_coords_path, "farm_temporal_data.qs")) %>%
    rename(farm_ID = farm_id)

  # Build farm_temporal_data: filter to production window and append Npop
  message("Building farm temporal data with population timeseries...")

  farm_temporal_data <- {
    ts_by_farm <- farm_temporal_data_raw %>% arrange(farm_ID) %>% split(.$farm_ID)
    st_by_farm <- farm_static_data       %>% arrange(farm_ID) %>% split(.$farm_ID)

    purrr::map2_dfr(st_by_farm, ts_by_farm, function(st, ts) {
      cbind(
        ts %>% filter(day >= st$t_start & day <= st$t_end),
        data.frame(Npop = generate_pop(
          harvest_n = st$harvest_n,
          mort      = all_params["mortmyt"],
          times     = c(t_start = st$t_start, t_end = st$t_end, dt = 1)
        ))
      )
    })
  }
  rm(farm_temporal_data_raw)
  qsave(farm_temporal_data, farm_temporal_data_file) 
}

farm_IDs <- farm_static_data %>% pull(farm_ID) %>% unique() %>% sort() %>% 
  sample(farm_bite_size)

# Sample of farms (10 per country) used for full daily timeseries outputs
# farm_IDs_sample <- farm_static_data %>%
#   group_by(country) %>%
#   slice_sample(n = 10) %>%
#   ungroup() %>%
#   pull(farm_ID) %>%
#   sort()
# qsave(farm_IDs_sample, file.path(out_root, "farm_IDs_sample.qs"))
farm_IDs_sample <- qread(file.path(out_root, "farm_IDs_sample.qs"))

message(sprintf(
  "Data loaded: %d farm IDs, %d feeds, %d sample farms.",
  length(farm_IDs), length(feed_names), length(farm_IDs_sample)
))


# ==========================================================================================================================================================
# 2. SECTION: farm_runs_full
#    Per-fish daily timeseries from Monte-Carlo runs, no population applied.
#    One file per farm_ID × feed_name.
# ==========================================================================================================================================================

message("\n--- SECTION: farm_runs_full ---")

# Build full grid of (farm_ID, feed_index) and filter to not-yet-done
all_farm_feed_combos <- expand.grid(
  farm_ID    = farm_IDs,
  feed_index = seq_along(feed_names),
  stringsAsFactors = FALSE
) %>%
  mutate(
    feed_name = feed_names[feed_index],
    fname     = file.path(
      dir_farm_runs_full, feed_name, 
      paste0("farm_runs_full_", fixnum(farm_ID), "_", feed_name, ".qs")
    ),
    fname_2   = file.path(
      bigdir_farm_runs_full, feed_name, 
      paste0("farm_runs_full_", fixnum(farm_ID), "_", feed_name, ".qs")
    ),
    fname_exists = file.exists(fname),
    fname_2_exists = file.exists(fname_2)
  ) %>%
  filter(!file.exists(fname) & !file.exists(fname_2))

message(sprintf(
  "%d / %d farm_runs_full files to generate - approx %.1f GB needed",
  nrow(all_farm_feed_combos),
  length(farm_IDs) * length(feed_names),
  1.2 * 0.085 * nrow(all_farm_feed_combos)
))

tic("farm_runs_full")

if (nrow(all_farm_feed_combos) != 0) {
  # plan(multisession, workers = parallelly::availableCores() - 1)
  plan(multisession, workers = 12)

  # Prepare static and temporal lookup lists
  static_by_farm  <- farm_static_data %>% filter(farm_ID %in% farm_IDs) %>% arrange(farm_ID) %>% split(.$farm_ID)
  temporal_by_farm <- farm_temporal_data %>% filter(farm_ID %in% farm_IDs) %>% arrange(farm_ID, day) %>% split(.$farm_ID)

  future_walk(
    seq_len(nrow(all_farm_feed_combos)),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE,
    function(i) {
      row     <- all_farm_feed_combos[i, ]
      fid     <- as.character(row$farm_ID)
      fi      <- row$feed_index
      outfile <- row$fname

      st <- static_by_farm[[fid]]
      ts <- temporal_by_farm[[fid]]

      result <- farm_growth_full(
        species_params = all_params,
        water_temp     = ts$temp_c,
        feed_params    = feed_params[[fi]],
        times          = c(t_start = st$t_start, t_end = st$t_end, dt = 1),
        output_vars    = c("weight", stats_measures),
        MC_pop         = inds_per_farm
      )

      qsave(result, outfile)
    }
  )

  plan(sequential)
  rm(static_by_farm, temporal_by_farm)
}

toc()
gc()


# ==========================================================================================================================================================
# 3. SECTION: fish_ends
#    Harvest-end per-fish summary. Loads farm_runs_full files within the map.
#    One file per farm_ID × feed_name.
# ==========================================================================================================================================================

message("\n--- SECTION: fish_ends ---")

fish_ends_combos <- expand.grid(
  farm_ID    = farm_IDs,
  feed_index = seq_along(feed_names),
  stringsAsFactors = FALSE
) %>%
  mutate(
    feed_name  = feed_names[feed_index],
    src_fname_1  = file.path(dir_farm_runs_full, feed_name, paste0("farm_runs_full_", fixnum(farm_ID), "_", feed_name, ".qs")),
    src_fname_2  = file.path(bigdir_farm_runs_full, feed_name, paste0("farm_runs_full_", fixnum(farm_ID), "_", feed_name, ".qs")),
    src_fname = case_when(
      file.exists(src_fname_1) ~ src_fname_1,
      file.exists(src_fname_2) ~ src_fname_2,
      TRUE ~ NA
    ),
    out_fname  = file.path(dir_fish_ends, feed_name, paste0("fish_ends_", fixnum(farm_ID), "_", feed_name, ".qs"))
  ) %>%
  filter(!file.exists(out_fname) & file.exists(src_fname))

message(sprintf(
  "%d / %d fish_ends files to generate.",
  nrow(fish_ends_combos),
  length(farm_IDs) * length(feed_names)
))

tic("fish_ends")

walk(
  seq_len(nrow(fish_ends_combos)),
  .progress = TRUE,
  function(i) {
    row       <- fish_ends_combos[i, ]
    farm_id   <- row$farm_ID
    feed_nm   <- row$feed_name
    src       <- row$src_fname
    outfile   <- row$out_fname

    ls_full <- qread(src)

    n_fish <- nrow(ls_full[["weight"]])
    n_days <- ncol(ls_full[["weight"]])

    harvest_weight <- ls_full[["weight"]][, n_days]
    measure_names  <- setdiff(names(ls_full), c("days", "weight"))

    sums <- lapply(measure_names, function(m) {
      matrixStats::rowSums2(ls_full[[m]], na.rm = TRUE)
    }) %>% setNames(measure_names)

    fe_weight <- data.frame(
      farm_ID       = as.integer(farm_id),
      feed          = as.factor(feed_nm),
      fish          = seq_len(n_fish),
      measure       = "weight",
      value         = harvest_weight,
      value_perbiom = 1
    )

    fe_outcomes <- purrr::map_dfr(measure_names, function(m) {
      data.frame(
        farm_ID       = as.integer(farm_id),
        feed          = as.factor(feed_nm),
        fish          = seq_len(n_fish),
        measure       = as.factor(m),
        value         = sums[[m]],
        value_perbiom = sums[[m]] / harvest_weight
      )
    })

    qsave(bind_rows(fe_weight, fe_outcomes), outfile)
  }
)

toc()
gc()


# ==========================================================================================================================================================
# 4. SECTION: farm_runs
#    Population-weighted daily timeseries, for sample farms only.
#    Loads farm_runs_full files within the map.
#    One file per farm_ID × feed_name.
# ==========================================================================================================================================================

message("\n--- SECTION: farm_runs ---")

farm_runs_combos <- expand.grid(
  farm_ID    = farm_IDs_sample,
  feed_index = seq_along(feed_names),
  stringsAsFactors = FALSE
) %>%
  mutate(
    feed_name = feed_names[feed_index],
    src_fname_1  = file.path(dir_farm_runs_full, feed_name, paste0("farm_runs_full_", fixnum(farm_ID), "_", feed_name, ".qs")),
    src_fname_2  = file.path(bigdir_farm_runs_full, feed_name, paste0("farm_runs_full_", fixnum(farm_ID), "_", feed_name, ".qs")),
    src_fname = case_when(
      file.exists(src_fname_1) ~ src_fname_1,
      file.exists(src_fname_2) ~ src_fname_2,
      TRUE ~ NA
    ),
    out_fname  = file.path(dir_farm_runs, feed_name, paste0("farm_runs_", fixnum(farm_ID), "_", feed_name, ".qs"))
  ) %>%
  filter(
    !file.exists(out_fname) &
      file.exists(src_fname)
  )

message(sprintf(
  "%d / %d farm_runs files to generate.",
  nrow(farm_runs_combos),
  length(farm_IDs_sample) * length(feed_names)
))

# Prepare Npop lookup for sample farms: daily Npop for each farm
ts_sample <- farm_temporal_data %>%
  filter(farm_ID %in% farm_IDs_sample) %>%
  arrange(farm_ID, day) %>%
  select(farm_ID, day, Npop) %>%
  split(.$farm_ID)

tic("farm_runs")

walk(
  seq_len(nrow(farm_runs_combos)),
  .progress = TRUE,
  function(i) {
    row       <- farm_runs_combos[i, ]
    farm_id   <- row$farm_ID
    feed_nm   <- row$feed_name
    src       <- row$src_fname
    outfile   <- row$out_fname

    df_decomposed <- qread(src)
    df_ts         <- ts_sample[[as.character(farm_id)]] %>%
      mutate(
        feed = as.factor(feed_nm),
        day  = day - min(day) + 1
      ) %>%
      select(farm_ID, feed, day, Npop)

    fa_pop <- decomposed_to_tidy(df_decomposed) %>%
      full_join(df_ts, by = c("prod_t" = "day")) %>%
      mutate(value = value * Npop) %>%
      select(-Npop)

    fa_weight <- fa_pop %>%
      filter(measure == "weight") %>%
      select(-measure) %>%
      rename(weight = value)

    fa_outcomes <- fa_pop %>%
      left_join(fa_weight, by = join_by(farm_ID, feed, fish, prod_t)) %>%
      mutate(value_perbiom = value / weight) %>%
      group_by(farm_ID, feed, prod_t, measure) %>%
      reframe(
        mean_value         = mean(value),
        sd_value           = sd(value),
        mean_value_perbiom = mean(value_perbiom),
        sd_per_biom        = sd(value_perbiom),
        n                  = as.integer(inds_per_farm)
      )

    qsave(fa_outcomes, outfile)
  }
)

toc()
rm(ts_sample)
gc()


# ==========================================================================================================================================================
# 5. SECTION: farm_ends
#    Population-weighted harvest summaries (all farms).
#    Loads fish_ends files within the map.
#    One file per farm_ID × feed_name.
# ==========================================================================================================================================================

message("\n--- SECTION: farm_ends ---")

farm_ends_combos <- expand.grid(
  farm_ID    = farm_IDs,
  feed_index = seq_along(feed_names),
  stringsAsFactors = FALSE
) %>%
  mutate(
    feed_name = feed_names[feed_index],
    src_fname  = file.path(dir_fish_ends, feed_name, paste0("fish_ends_", fixnum(farm_ID), "_", feed_name, ".qs")),,
    out_fname = file.path(
      dir_farm_ends, feed_name,
      paste0("farm_ends_", fixnum(farm_ID), "_", feed_name, ".qs")
    )
  ) %>%
  filter(!file.exists(out_fname) & file.exists(src_fname))

message(sprintf(
  "%d / %d farm_ends files to generate.",
  nrow(farm_ends_combos),
  length(farm_IDs) * length(feed_names)
))

# Pre-build harvest-day Npop lookup: one row per farm (last day = harvest)
harvest_npop <- farm_temporal_data %>% 
  group_by(farm_ID) %>%
  slice_max(day, n = 1) %>%
  ungroup() %>%
  select(farm_ID, Npop)

tic("farm_ends")

walk(
  seq_len(nrow(farm_ends_combos)),
  .progress = TRUE,
  function(i) {
    row     <- farm_ends_combos[i, ]
    farm_id <- row$farm_ID
    feed_nm <- row$feed_name
    src     <- row$src_fname
    outfile <- row$out_fname

    df_fish_end <- qread(src)
    npop        <- harvest_npop$Npop[harvest_npop$farm_ID == farm_id]

    fa_pop <- df_fish_end %>%
      mutate(value = value * npop) %>%
      select(-value_perbiom)

    fa_weight <- fa_pop %>%
      filter(measure == "weight") %>%
      select(-measure) %>%
      rename(weight = value)

    fa_outcomes <- fa_pop %>%
      left_join(fa_weight, by = join_by(farm_ID, feed, fish)) %>%
      mutate(value_perbiom = value / weight) %>%
      group_by(farm_ID, feed, measure) %>%
      reframe(
        mean_value         = mean(value),
        sd_value           = sd(value),
        mean_value_perbiom = mean(value_perbiom),
        sd_per_biom        = sd(value_perbiom),
        n                  = as.integer(inds_per_farm)
      )

    qsave(fa_outcomes, outfile)
  }
)

toc()

rm(harvest_npop)
gc()

# Move large files
message("\nMoving farm_runs_full files to big storage...")
all_farm_feed_combos <- expand.grid(
  farm_ID    = farm_IDs,
  feed_index = seq_along(feed_names),
  stringsAsFactors = FALSE
) %>%
  mutate(
    feed_name = feed_names[feed_index],
    fname     = file.path(
      dir_farm_runs_full, feed_name, 
      paste0("farm_runs_full_", fixnum(farm_ID), "_", feed_name, ".qs")
    ),
    fname_2   = file.path(
      bigdir_farm_runs_full, feed_name, 
      paste0("farm_runs_full_", fixnum(farm_ID), "_", feed_name, ".qs")
    ),
    fname_exists = file.exists(fname),
    fname_2_exists = file.exists(fname_2)
  ) %>% 
  filter(fname_exists)

file.rename(
  all_farm_feed_combos$fname, 
  all_farm_feed_combos$fname_2
)

message("\nAll pipeline sections complete.\n")
