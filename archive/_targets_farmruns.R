suppressPackageStartupMessages(suppressWarnings({
  library(targets)
  library(crew)
  library(qs)
  library(tidyverse)
  library(lme4)
}))

tar_option_set(
  format = "qs", 
  controller = crew::crew_controller_local(
    workers = 4, 
    seconds_idle = 120,
    crashes_max = 10
    ),
  error = "abridge",
  workspace_on_error = TRUE,
  garbage_collection = 2
)

tar_source(
  files = c(
    here::here("src/dirs.R"),
    here::here("src/functions.R"),
    here::here("src/model_functions_new.R")
  )
)

# Dirs
output_path <- here::here("outputs")
output_farm_data_path <- file.path(output_path, "farm_data")

inds_per_farm <- 2000

list(
# Get previously saved data ------------------------------------------------------------------------------------------------------------------------
  tar_target(species_params_file, file.path(output_species_data_path, "species_params.qs"), format = "file"),
  tar_target(pop_params_file, file.path(output_species_data_path, "pop_params.qs"), format = "file"),
  tar_target(farm_temporal_data_file, file.path(input_farm_coords_path, "farm_temporal_data.qs"), format = "file"),
  tar_target(farm_static_data_file, file.path(input_farm_coords_path, "farm_static_data_w_stocking.qs"), format = "file"),

  tar_target(feed_params_file_MD, file.path(output_species_data_path, "feed_params_MD.qs"), format = "file"),
  tar_target(feed_params_file_PD, file.path(output_species_data_path, "feed_params_PD.qs"), format = "file"),
  tar_target(feed_params_file_AD, file.path(output_species_data_path, "feed_params_AD.qs"), format = "file"),

  tar_target(all_params, c(qread(species_params_file), qread(pop_params_file))),

  tar_target(
    farm_IDs,
    description = {"List all the unique farm_IDs being run"},
    command = {
      farm_static_data_file %>% 
        qread() %>% 
        pull(farm_id) %>% 
        unique() %>% 
        sort() #%>% head(500)
    }
  ),

  tar_target(
    farm_IDs_chunked,
    description = {"Chunk farm_IDs to reduce the number of total objects"},
    command = {split(farm_IDs, ceiling(seq_along(farm_IDs) / 10))}
  ),

  tar_target(
    farm_IDs_sample,
    description = {"Subset a sample of farms to keep the full timeseries of for supplementary plotting"},
    command = {
      farm_static_data_file %>% 
        qread() %>% 
        sf::st_drop_geometry() %>% 
        group_by(country) %>% 
        slice_sample(n = 5) %>% 
        ungroup() %>% 
        pull(farm_id) %>% 
        sort()
    }
  ),

  tar_target(
    feed_params, 
    description = {"Specifics of all the feeds being run"},
    command = {c(qread(feed_params_file_MD), qread(feed_params_file_PD), qread(feed_params_file_AD))}
  ), 
  tar_target(feed_names, names(feed_params)), 

  tar_target(
    farm_static_data, 
    description = {"Location-specific data for all farms"},
    command = {
      qread(farm_static_data_file) %>% 
        sf::st_drop_geometry() %>% 
        rename(farm_ID = farm_id)
  }),

  tar_target(
    farm_temporal_data,
    description = {"Timeseries data for all farms"},
    command = {
      ts_data <- qread(farm_temporal_data_file)%>% 
        rename(farm_ID = farm_id) %>% 
        arrange(farm_ID) %>% 
        split(.$farm_ID)

      st_data <- farm_static_data %>% 
        arrange(farm_ID) %>% 
        split(.$farm_ID)

      purrr::map2_dfr(
        st_data, 
        ts_data, 
        function(st, ts) {
          cbind(
            ts %>% 
              filter(day >= st$t_start & day <= st$t_end),
            data.frame(Npop = generate_pop(
                harvest_n = st$harvest_n,
                mort = all_params['mortmyt'],
                times = c(t_start = st$t_start, t_end = st$t_end, dt = 1)
            ))
          )
        })
    }),

  tar_target(
    stats_measures, 
    description = {"Measures that the farm_runs will output"},
    c("food_prov", "total_excr", "total_uneat", "nitrogen_excr", "nitrogen_uneat", "total_nitrogen", "carbon_excr", "carbon_uneat", "total_carbon")
  ),

  tar_target(
    farm_runs_full,
    description = {"Basic farm runs, with the Monte-Carlo sampling but without applying individual farm populations. Seperate fish are preserved."},
    command = {
      static <- farm_static_data %>% 
        filter(farm_ID %in% farm_IDs_chunked[[1]]) %>% 
        arrange(farm_ID) %>% 
        split(.$farm_ID)
      temporal <- farm_temporal_data %>% 
        filter(farm_ID %in% farm_IDs_chunked[[1]]) %>% 
        arrange(farm_ID, day) %>% 
        split(.$farm_ID)

      # st <- static[[1]]
      # ts <- temporal[[1]]
      map2(static, temporal, function(st, ts) {
        farm_growth_full(
          species_params = all_params,
          water_temp = ts$temp_c,
          feed_params = feed_params[[1]],
          times = c(t_start = st$t_start, t_end = st$t_end, dt = 1),
          output_vars = c("weight", stats_measures),
          MC_pop = inds_per_farm
        )
      }) %>% 
        setNames(farm_IDs_chunked[[1]])
    },
    pattern = cross(farm_IDs_chunked, map(feed_params, feed_names)),
    iteration = "list"
  ),

  tar_target(
    fish_ends,
    description = {"Harvest values only. Farm populations are not applied, and seperate fish are preserved."},
    command = {
      # ls_full <- farm_runs_full[[1]]
      # farm_id <- farm_IDs_chunked[[1]][[1]]
      map2_dfr(farm_runs_full, farm_IDs_chunked[[1]], function(ls_full, farm_id) {
        n_fish <- nrow(ls_full[["weight"]])
        n_days <- ncol(ls_full[["weight"]])

        # Harvest weight of each fish (last time-step column)
        harvest_weight <- ls_full[["weight"]][, n_days]

        # Measure names to sum (exclude 'days' and 'weight')
        measure_names <- setdiff(names(ls_full), c("days", "weight"))

        # Summed outputs over production period for each fish
        sums <- lapply(measure_names, function(m) matrixStats::rowSums2(ls_full[[m]], na.rm = T)) %>%
          setNames(measure_names)

        # Harvest weight row (value_perbiom = 1 by definition)
        fe_weight <- data.frame(
          farm_ID       = as.integer(farm_id),
          feed          = as.factor(feed_names),
          fish          = seq_len(n_fish),
          measure       = "weight",
          value         = harvest_weight,
          value_perbiom = 1
        )

        # Summed outcomes per fish, normalised by harvest biomass
        fe_outcomes <- map_dfr(measure_names, function(m) {
          data.frame(
            farm_ID       = as.integer(farm_id),
            feed          = as.factor(feed_names),
            fish          = seq_len(n_fish),
            measure       = m,
            value         = sums[[m]],
            value_perbiom = sums[[m]] / harvest_weight
          )
        })
        bind_rows(fe_weight, fe_outcomes)
      })
    },
    pattern = map(farm_runs_full, cross(farm_IDs_chunked, feed_names))
  ),

  tar_target(
    farm_runs,
    description = {"All daily values from a small sample of farms. Farm populations are applied and fish are combined into farm mean and sd."},
    command = {
      ls_decomposed <- farm_runs_full[names(farm_runs_full) %in% farm_IDs_sample]
      ls_ts <- farm_temporal_data %>% 
        filter(farm_ID %in% names(ls_decomposed)) %>% 
        mutate(feed = feed_names) %>% 
        arrange(farm_ID, day) %>% 
        select(farm_ID, feed, day, Npop) %>% 
        arrange(farm_ID) %>% 
        split(.$farm_ID)

      # df_decomposed <- ls_decomposed[[1]]
      # df_ts <- ls_ts[[1]]
      if (length(ls_decomposed) > 0) {
        map2_dfr(ls_decomposed, ls_ts, function(df_decomposed, df_ts) {
          fa_pop <- decomposed_to_tidy(df_decomposed) %>% 
            full_join(
              df_ts %>% 
                mutate(day = day - min(day) + 1),
              by = c("prod_t" = "day")
            ) %>% 
            mutate(value = value * Npop) %>% 
            select(-Npop)
          
          # Harvest biomass
          fa_weight <- fa_pop %>% 
            filter(measure == "weight") %>% 
            select(-measure) %>% 
            rename(weight = value)

          # Summed outputs over production period for each fish (with population)
          fa_outcomes <- fa_pop %>%
            left_join(fa_weight, by = join_by(farm_ID, feed, fish, prod_t)) %>% 
            mutate(value_perbiom = value/weight) %>% 
            group_by(farm_ID, feed, prod_t, measure) %>% 
            reframe(
              mean_value = mean(value),
              sd_value = sd(value),
              mean_value_perbiom = mean(value_perbiom),
              sd_per_biom = sd(value_perbiom),
              n = as.integer(inds_per_farm)
            )
        })
      }
    },
    pattern = map(farm_runs_full, cross(farm_IDs_chunked, feed_names))
  ),

  tar_target(
    farm_ends,
    description = {"Harvest values only. Farm populations are applied and fish are combined into farm mean and sd."},
    command = {
      ls_fish_ends <- fish_ends %>% 
        arrange(farm_ID) %>% 
        split(.$farm_ID)
      ls_ts <- farm_temporal_data %>% 
        filter(farm_ID %in% names(ls_fish_ends)) %>% 
        group_by(farm_ID) %>% 
        slice_max(day) %>% 
        select(farm_ID, Npop) %>% 
        arrange(farm_ID) %>% 
        split(.$farm_ID)

      # df_fish_end <- ls_fish_ends[[1]]
      # df_ts <- ls_ts[[1]]
      map2_dfr(ls_fish_ends, ls_ts, function(df_fish_end, df_ts) {
        fa_pop <- df_fish_end %>% 
          mutate(value = value * df_ts$Npop) %>% 
          select(-value_perbiom)

        # Harvest biomass
        fa_weight <- fa_pop %>% 
          filter(measure == "weight") %>% 
          select(-measure) %>% 
          rename(weight = value)

        # Summed outputs over production period for each fish (with population)
        fa_outcomes <- fa_pop %>%
          left_join(fa_weight, by = join_by(farm_ID, feed, fish)) %>% 
          mutate(value_perbiom = value/weight) %>% 
          group_by(farm_ID, feed, measure) %>% 
          reframe(
            mean_value = mean(value),
            sd_value = sd(value),
            mean_value_perbiom = mean(value_perbiom),
            sd_per_biom = sd(value_perbiom),
            n = as.integer(inds_per_farm)
          )
      })
    },
    pattern = fish_ends
  )
)
