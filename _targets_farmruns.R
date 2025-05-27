library(targets)
library(crew)

tar_option_set(
  packages = c("stringr", "magrittr", "tidyr", "arrow", "dplyr", "tibble"), 
  format = "qs", 
  controller = crew_controller_local(workers = 16),
  workspace_on_error = TRUE
)

# Helper function to process each matrix
farm_to_cohort <- function(matrix, time_offset = 0) {
  matrix %>% 
    as.data.frame() %>% 
    rename(t = V1, mean = V2, sd = V3) %>%
    mutate(t = t + time_offset)
}

# Process each time period (current, +365 days, +730 days)
combine_cohorts <- function(lst) {
  bind_rows(farm_to_cohort(lst[[i]], 0),
            farm_to_cohort(lst[[i]], 365),
            farm_to_cohort(lst[[i]], 730))
}

list(
# Prepare parameters ----------------------------------------------------------------------------------------------
  tar_target(
    tar_farm_IDs,
    farm_ts_data %>% distinct(farm_ID) %>% pull(farm_ID) #%>% sample(500)
  ),

  tar_target(
    tar_common_params, 
    command = c(
      't_start' = farm_coords$t_start[farm_coords$farm_ID == tar_farm_IDs],
      't_end' = farm_coords$t_end[farm_coords$farm_ID == tar_farm_IDs],
      'nruns' = 5000
    ),
    pattern = tar_farm_IDs
  ),
  
  tar_target(
    tar_farm_ts,
    command = {
      farm_ts_data %>%
        dplyr::filter(farm_ID == tar_farm_IDs) %>%
        dplyr::filter(between(day, tar_common_params['t_start'], tar_common_params['t_end']))
    },
    pattern = map(tar_farm_IDs, tar_common_params)
  ),
  
  tar_target(
    tar_static_data, 
    file.path(input_farm_coords_path, str_c(this_species, "_locations_w_temps.qs")) %>% 
      qs::qread() %>% 
      distinct(farm_id, tonnes_per_farm) %>% 
      mutate(tonnes_per_farm = tonnes_per_farm %>% units::set_units("t") %>% units::set_units("g") %>% units::drop_units())
  ),
  tar_target(
    tar_harvest_size,
    farm_harvest_size$weight[farm_harvest_size$farm_ID == tar_farm_IDs],
    pattern = tar_farm_IDs
  ),
  tar_target(
    tar_N_pop, 
    command = {
      generate_pop(
        harvest_n = (tar_static_data$tonnes_per_farm[tar_static_data$farm_id == tar_farm_IDs]/tar_harvest_size),
        mort = pop_params['mortmyt'],
        times = c(tar_common_params['t_start'], tar_common_params['t_end'], 'dt' = 1)
      )
    }, 
    pattern = map(tar_farm_IDs, tar_harvest_size, tar_common_params)
  ),
  
# Run growth ------------------------------------------------------------------------------------------------------
  tar_target(
    tar_farmrun_reference, 
    command = {
      farm_growth(
        pop_params = pop_params,
        species_params = species_params,
        water_temp = tar_farm_ts$temp_c,
        feed_params = feed_params[['reference']],
        times = c(tar_common_params['t_start'], tar_common_params['t_end'], 'dt' = 1),
        N_pop = unname(tar_N_pop),
        nruns = tar_common_params['nruns']
      )
    },
    pattern = map(tar_farm_ts, tar_farm_IDs, tar_N_pop, tar_common_params)
  ),
  tar_target(
    tar_farmrun_past, 
    command = {
      farm_growth(
        pop_params = pop_params,
        species_params = species_params,
        water_temp = tar_farm_ts$temp_c,
        feed_params = feed_params[['past']],
        times = c(tar_common_params['t_start'], tar_common_params['t_end'], 'dt' = 1),
        N_pop = tar_N_pop,
        nruns = tar_common_params['nruns']
      )
    },
    pattern = map(tar_farm_ts, tar_farm_IDs, tar_N_pop, tar_common_params)
  ),
  tar_target(
    tar_farmrun_future, 
    command = {
      farm_growth(
        pop_params = pop_params,
        species_params = species_params,
        water_temp = tar_farm_ts$temp_c,
        feed_params = feed_params[['future']],
        times = c(tar_common_params['t_start'], tar_common_params['t_end'], 'dt' = 1),
        N_pop = tar_N_pop,
        nruns = tar_common_params['nruns']
      )
    },
    pattern = map(tar_farm_ts, tar_farm_IDs, tar_N_pop, tar_common_params)
  ),
  
  tar_target(
    tar_farmrun_comparisons,
    command = {
      # Convert farms to cohorts (reference feed)
      ref1 <- lapply(tar_farmrun_reference, farm_to_cohort, time_offset = 0)
      ref2 <- lapply(tar_farmrun_reference, farm_to_cohort, time_offset = 365)
      ref3 <- lapply(tar_farmrun_reference, farm_to_cohort, time_offset = 730)
      
      # Get the t values of cohort 2
      lims <- unique(ref2[[1]]$t)
      
      ref <- purrr::map(1:length(stat_names), function(st) {
        bind_rows(ref1[[st]], ref2[[st]], ref3[[st]]) %>% 
          mutate(sd = sd/mean, 
                 feed = factor("reference", levels = c("reference", "past", "future"))) %>%
          group_by(feed, t) %>%
          reframe(mean = sum(mean),
                  sd = sum(sd)*mean) %>% 
          filter(t %in% lims)
      })
      
      # Convert farms to cohorts (future feed)
      ref1 <- lapply(tar_farmrun_future, farm_to_cohort, time_offset = 0)
      ref2 <- lapply(tar_farmrun_future, farm_to_cohort, time_offset = 365)
      ref3 <- lapply(tar_farmrun_future, farm_to_cohort, time_offset = 730)
      
      # Get the t values of cohort 2
      lims <- unique(ref2[[1]]$t)
      
      fut <- purrr::map(1:length(stat_names), function(st) {
        bind_rows(ref1[[st]], ref2[[st]], ref3[[st]]) %>% 
          mutate(sd = sd/mean, 
                 feed = factor("future", levels = c("reference", "past", "future"))) %>%
          group_by(feed, t) %>%
          reframe(mean = sum(mean),
                  sd = sum(sd)*mean) %>% 
          filter(t %in% lims)
      })
      
      # Convert farms to cohorts (past feed)
      ref1 <- lapply(tar_farmrun_past, farm_to_cohort, time_offset = 0)
      ref2 <- lapply(tar_farmrun_past, farm_to_cohort, time_offset = 365)
      ref3 <- lapply(tar_farmrun_past, farm_to_cohort, time_offset = 730)
      
      # Get the t values of cohort 2
      lims <- unique(ref2[[1]]$t)
      
      pas <- purrr::map(1:length(stat_names), function(st) {
        bind_rows(ref1[[st]], ref2[[st]], ref3[[st]]) %>% 
          mutate(sd = sd/mean, 
                 feed = factor("past", levels = c("reference", "past", "future"))) %>%
          group_by(feed, t) %>%
          reframe(mean = sum(mean),
                  sd = sum(sd)*mean) %>% 
          filter(t %in% lims)
      })
      
      # Put them all together
      purrr::map(1:length(stat_names), function(st) {
        bind_rows(ref[[st]], pas[[st]], fut[[st]]) %>% 
          mutate(farm_ID = tar_farm_IDs)
      })
    },
    pattern = map(tar_farm_IDs, tar_farmrun_reference, tar_farmrun_past, tar_farmrun_future)
  )
)
