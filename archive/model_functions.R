### Main model functions here are modified from Baldan et al 2018 R package for aquaculture. 
### https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195732
### https://github.com/cran/RAC/tree/master/R

suppressPackageStartupMessages(suppressWarnings({
  library(qs)
  library(terra)
  library(readxl)
  library(matrixStats)
  library(dplyr)
  library(msm)
  library(reshape2)
  library(purrr)
  library(future)
  library(furrr)
}))

# Parameters definitions
# species_params['alpha']         [-] Feeding catabolism coefficient
# species_params['epsprot']       [J/gprot] Energy content of protein
# species_params['epslip']        [J/glip] Energy content of lipid
# species_params['epscarb']       [J/gcarb] Energy content of carbohydrate
# species_params['epsO2']         [J/gO2] Energy consumed by the respiration of 1g of oxygen
# species_params['pk']            [1/day] Temperature coefficient for the fasting catabolism
# species_params['k0']            [1/Celsius degree] Fasting catabolism at 0 Celsius degree
# species_params['m']             [-] Weight exponent for the anabolism
# species_params['n']             [-] Weight exponent for the catabolism
# species_params['betac']         [-] Shape coefficient for the H(Tw) function
# species_params['Tma']           [Celsius degree] Maximum feeding temperature  
# species_params['Toa']           [Celsius degree] Optimal feeding temperature
# species_params['Taa']           [Celsius degree] Lowest feeding temperature
# species_params['omega']         [gO2/g] Oxygen consumption - weight loss ratio
# species_params['a']             [J/gtissue] Energy content of fish tissue
# species_params['k']             [-] Weight exponent for energy content
# species_params['eff']           [-] Food ingestion efficiency
# species_params['fcr']           [-] Food conversion ratio

# Function to get N (g) from protein (g)
get_nitrogen <- function(P) {unname(P/6.25)}

# Function to get C (g) from protein, lipid and carbs (g)
get_carbon <- function(P, L, C) {
  carbon_1 <- L * 0.75
  carbon_2 <- C * 0.41
  mol_N_P  <- get_nitrogen(P=P)/14.007
  mol_C_P  <- mol_N_P * 3.7
  carbon_3 <- mol_C_P * 12.011
  unname(carbon_1 + carbon_2 + carbon_3)
}

# Function to generate a population timeseries from harvest population and fixed mortality rate
generate_pop <- function(harvest_n, mort, times) {
  
  ts <- seq(times['t_start'], times['t_end'], by = times['dt'])   # Integration times
  
  # Initial condition and vectors initialization
  N_pop <- rep(0, length(ts))                              # Initialize vector N_pop
  N_pop[1] <- harvest_n                                    # Impose harvest condition
  
  # for cycle that solves population ODE with Euler method
  for (t in 2:length(ts)){
    dN <- unname(mort*N_pop[t-1])                           # Individuals increment
    N_pop[t] <- N_pop[t-1]+dN#*times['dt']                  # Individuals at time t+1
    
    # # Taking out the management alterations for now
    # for (i in 1:length(manag[,1])) {  # For cycle that adjusts N_pop according with management strategies
    #   if (t==manag[i,1]) {              # if statement to check if it is the time to adjust N_pop
    #     N_pop[t+1]=N_pop[t]+manag[i,2]
    #   } 
    # } 
  }
  return(rev(N_pop))
}

# Function to generate relative feeding rate (temperature dependent)
feeding_rate <- function(water_temp, species_params) {
  exp(species_params['betac'] * (water_temp - species_params['Toa'])) * 
    ((species_params['Tma'] - water_temp)/(species_params['Tma'] - species_params['Toa']))^
    (species_params['betac'] * (species_params['Tma'] - species_params['Toa']))
}

# Function to calculate how much food will be provided based on need
food_prov_rate <- function(pop_params, water_temp, ing_pot, ing_pot_min, species_params) {
  # Use ifelse vectorization instead of individual if statements
  if (water_temp > species_params['Taa'] & water_temp < species_params['Tma']) {
    ing_pot * (1 + rnorm(1, pop_params['overFmean'], pop_params['overFdelta']))
  } else {
    ing_pot_min
   } # old formula: 0.25 * 0.066 * weight^0.75
}

# Function to apportion ingested feed into relevant components (uneaten feed, excreted faeces, assimilated feed)
app_feed <- function(provided, ingested, prop, macro, digestibility) {
  # Pre-compute common values and use vectorized operations
  provided_g <- provided * prop * macro
  ingested_g <- ingested * prop * macro
  assimilated_g <- ingested_g * digestibility
  
  # Return only necessary values in a numeric vector
  c(
    # provided = sumna(provided_g),
    # ingested = sumna(ingested_g),
    uneaten = sumna(provided_g - ingested_g),
    assimilated = sumna(assimilated_g),
    excreted = sumna(ingested_g - assimilated_g)
  )
}

# Main fish growth function
fish_growth <- function(pop_params, species_params, water_temp, feed_params, times, init_weight, ingmax) {
  # Pre-calculate array sizes
  n_days <- length(times['t_start']:times['t_end'])
  
  # Preallocate all vectors at once
  result <- matrix(0, nrow = n_days, ncol = 22)
  colnames(result) <- c('days', 'weight', 'dw', 'water_temp', 'T_response', 'P_excr', 
                        'L_excr', 'C_excr', 'P_uneat', 'L_uneat', 'C_uneat', 'food_prov', 
                        'food_enc', 'rel_feeding', 'ing_pot', 'ing_act', 'E_assim', 
                        'E_somat', 'anab', 'catab', 'O2', 'NH4')
  
  # Initialize first values
  result[, 'days'] <- (times['t_start']:times['t_end'])*times['dt']
  result[1, 'weight'] <- init_weight
  result[, 'water_temp'] <- water_temp
  
  # Main calculation loop
  for (i in 1:(n_days-1)) {
    # Temperature response and ingestion calculations
    result[i, 'rel_feeding'] <- feeding_rate(result[i, 'water_temp'], species_params)
    result[i, 'ing_pot'] <- ingmax * (result[i, 'weight']^species_params['m']) * result[i, 'rel_feeding']
    
    # Food provision and ingestion
    result[i, 'food_prov'] <- food_prov_rate(
      pop_params = pop_params, 
      water_temp = result[i, 'water_temp'],
      ing_pot = result[i, 'ing_pot'],
      ing_pot_min = ingmax * (result[i, 'weight']^species_params['m']) * feeding_rate(species_params['Taa'], species_params),
      species_params
    )
    result[i, 'food_enc'] <- species_params['eff'] * result[i, 'food_prov']
    result[i, 'ing_act'] <- min(result[i, 'food_enc'], result[i, 'ing_pot'])
    
    # Energy calculations
    result[i, 'E_somat'] <- species_params['a'] * result[i, 'weight']^species_params['k']
    
    # Process feed components - vectorized operations
    app_carbs <- app_feed(
      provided = result[i, 'food_prov'], ingested = result[i, 'ing_act'], 
      prop = feed_params[['Carbohydrates']]$proportion, 
      macro = feed_params[['Carbohydrates']]$macro, 
      digestibility = feed_params[['Carbohydrates']]$digest
    )
    app_lipids <- app_feed(
      result[i, 'food_prov'], result[i, 'ing_act'],
      feed_params[['Lipids']]$proportion,
      feed_params[['Lipids']]$macro,
      feed_params[['Lipids']]$digest
    )
    app_proteins <- app_feed(
      result[i, 'food_prov'], result[i, 'ing_act'],
      feed_params[['Proteins']]$proportion,
      feed_params[['Proteins']]$macro,
      feed_params[['Proteins']]$digest
    )
    
    # Excretion and waste values
    result[i, c('C_excr', 'L_excr', 'P_excr')] <- c(app_carbs['excreted'], app_lipids['excreted'], app_proteins['excreted'])
    result[i, c('C_uneat', 'L_uneat', 'P_uneat')] <- c(app_carbs['uneaten'], app_lipids['uneaten'], app_proteins['uneaten'])

    # Energy assimilation
    result[i, 'E_assim'] <- 
      app_carbs['assimilated'] * species_params['epscarb'] +
      app_lipids['assimilated'] * species_params['epslip'] +
      app_proteins['assimilated'] * species_params['epsprot']
    
    # Temperature response and metabolism
    result[i, 'T_response'] <- exp(species_params['pk'] * result[i, 'water_temp'])
    result[i, 'anab'] <- result[i, 'E_assim'] * (1 - species_params['alpha'])
    result[i, 'catab'] <- species_params['epsO2'] * species_params['k0'] * 
      result[i, 'T_response'] * (result[i, 'weight']^species_params['n']) * 
      species_params['omega']
    
    # O2 and NH4 calculations
    result[i, 'O2'] <- result[i, 'catab'] / species_params['epsO2']
    result[i, 'NH4'] <- result[i, 'O2'] * 0.06
    
    # Weight calculations
    result[i, 'dw'] <- (result[i, 'anab'] - result[i, 'catab']) / result[i, 'E_somat']
    result[i + 1, 'weight'] <- result[i, 'weight'] + result[i, 'dw'] * times['dt']
  } 
  result
}

# Farm growth function - applies the fish growth function over a Monte-Carlo sampled population
farm_growth <- function(pop_params, species_params, feed_params, water_temp, times, N_pop, nruns){
    
  days <- (times['t_start']:times['t_end'])*times['dt']
  
  # Generate all random values upfront
  init_weights <- rnorm(nruns, mean = pop_params['meanW'], sd = pop_params['deltaW'])
  ingmaxes <- rnorm(nruns, mean = pop_params['meanImax'], sd = pop_params['deltaImax'])
  
  # Run parallel simulation for individuals
  mc_results <- future_map2(init_weights, ingmaxes, function(init_w, ing_m) {
    mat <- fish_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = water_temp,
      feed_params = feed_params,
      times = times,
      init_weight = init_w,
      ingmax = ing_m
    ) %>% unname()
  })
  
  stat_names <- c("days", "weight", "dw", "water_temp", "T_response", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "food_prov", "food_enc", "rel_feeding", "ing_pot", "ing_act", "E_assim", "E_somat", "anab", "catab", "O2", "NH4")

  # Consolidate all individuals into a farm (with population = nruns) - individuals are rows, timesteps are columns
  all_results <- lapply(1:length(stat_names), function(col_idx) {
    t(
      sapply(mc_results, function(mat_idx) {
        mat_idx[, col_idx]
      })
    )
  }) %>% 
    setNames(stat_names)
  
  # Some stats need to be summed/added
  all_results[["total_excr"]] <- all_results[["P_excr"]] + all_results[["L_excr"]] + all_results[["C_excr"]]
  all_results[["total_uneat"]] <- all_results[["P_uneat"]] + all_results[["L_uneat"]] + all_results[["C_uneat"]]
  all_results[["metab"]] <- all_results[["anab"]] - all_results[["catab"]]
  all_results[["biomass"]] <- all_results[["weight"]] # weight is individual weight, biomass will be farm biomass
  
  # Average across individuals to get mean and sd for the whole farm
  all_results <- lapply(2:length(all_results), function(col_idx) {
    cbind(
      colMeans(all_results[[col_idx]]), 
      matrixStats::colSds(all_results[[col_idx]])
    ) %>% 
      as.matrix() %>% unname()
  }) %>% 
    setNames(names(all_results)[2:length(names(all_results))])

  # Some stats need to be multiplied by the farm population (Npop)
  pop_names <- c("biomass", "dw", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "ing_act", "total_excr", "total_uneat", "O2", "NH4", "food_prov")
  for (stat_nm in pop_names) {
    all_results[[stat_nm]][,2] <- (all_results[[stat_nm]][,2]/all_results[[stat_nm]][,1])
    all_results[[stat_nm]][,1] <- all_results[[stat_nm]][,1] * N_pop[1:length(days)]
    all_results[[stat_nm]][,2] <- all_results[[stat_nm]][,1] * all_results[[stat_nm]][,2]
  }
  
  out_list <- lapply(1:length(all_results), function(col_idx) {
    cbind(days, all_results[[col_idx]]) %>% 
      as.matrix() %>% unname() 
  }) %>% setNames(paste0(names(all_results), "_stat"))
  
  return(out_list)
}

# This is identical to the farm_growth function (and still multiplies by population) except without the Monte-Carlo sampling of initial weights (all fish are uniform)
uni_farm_growth <- function(pop_params, species_params, feed_params, water_temp, times, N_pop, nruns = 1){
  
  new_pop_params <- pop_params
  new_pop_params['deltaW'] <- 0
  new_pop_params['deltaImax'] <- 0

  farm_growth(
    pop_params = new_pop_params, 
    species_params = species_params, 
    feed_params = feed_params, 
    water_temp = water_temp, 
    times = times, 
    N_pop = N_pop, 
    nruns = 1
  )
}

# This is identical to the above function BUT it does not condense each farm into a mean - it keeps the individual fish seperate
farm_growth_decomposed <- function(pop_params, species_params, feed_params, water_temp, times, care_about, nruns){
    
  days <- (times['t_start']:times['t_end'])*times['dt']

  # Generate all random values upfront
  init_weights <- rnorm(nruns, mean = pop_params['meanW'], sd = pop_params['deltaW'])
  ingmaxes <- rnorm(nruns, mean = pop_params['meanImax'], sd = pop_params['deltaImax'])
  
  # Run parallel simulation for individuals
  mc_results <- purrr::map2(init_weights, ingmaxes, function(init_w, ing_m) {
    mat <- fish_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = water_temp,
      feed_params = feed_params,
      times = times,
      init_weight = init_w,
      ingmax = ing_m
    ) %>% unname()
  })

  stat_names <- c("days", "weight", "dw", "water_temp", "T_response", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "food_prov", "food_enc", "rel_feeding", "ing_pot", "ing_act", "E_assim", "E_somat", "anab", "catab", "O2", "NH4")
  stat_names_2 <- c(stat_names, "total_excr", "total_uneat", "Nitrogen_excr", "Nitrogen_uneat", "Carbon_excr", "Carbon_uneat")

  # Consolidate all individuals into a farm (with population = nruns) - individuals are rows, timesteps are columns
  all_results <- map(1:length(stat_names), function(col_idx) {
    t(sapply(mc_results, function(mat_idx) {mat_idx[, col_idx]}))
    }) %>% 
    setNames(stat_names)
  
  all_results[["total_excr"]] <- all_results[["P_excr"]] + all_results[["L_excr"]] + all_results[["C_excr"]]
  all_results[["total_uneat"]] <- all_results[["P_uneat"]] + all_results[["L_uneat"]] + all_results[["C_uneat"]]
  all_results[["Nitrogen_excr"]] <- get_nitrogen(all_results[["P_excr"]])
  all_results[["Nitrogen_uneat"]] <- get_nitrogen(all_results[["P_uneat"]])
  all_results[["Carbon_excr"]] <- get_carbon(P = all_results[["P_excr"]], L = all_results[["L_excr"]], C = all_results[["C_excr"]])
  all_results[["Carbon_uneat"]] <- get_carbon(P = all_results[["P_uneat"]], L = all_results[["L_uneat"]], C = all_results[["C_uneat"]])
  all_results[["total_Carbon"]] <- all_results[["Carbon_excr"]] + all_results[["Carbon_uneat"]]
  all_results[["total_Nitrogen"]] <- all_results[["Nitrogen_excr"]] + all_results[["Nitrogen_uneat"]]

  # # Narrow down to the ones I actually care about
  all_results <- all_results[care_about]

  all_results <- map_dfr(1:length(all_results), function(col_idx) {
    res <- reshape2::melt(all_results[[col_idx]]) %>% mutate(measure = as.factor(care_about[col_idx]))
    colnames(res) <- c("fish", "prod_t", "value", "measure")
    res
  })
  return(all_results)
}

convert_feed_to_levels <- function(df) {
  df %>% 
    mutate(
      feed_level = case_when(str_detect(feed, "_min") ~ "min", str_detect(feed, "_max") ~ "max", T ~ "mean") %>% as.factor(),
      feed = str_remove_all(feed, "_min|_max") %>% as.factor()
    )
}
