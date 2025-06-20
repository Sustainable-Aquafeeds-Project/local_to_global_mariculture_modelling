# nolint start

# These packages are called here so that renv doesn't clean them
library(devtools)
library(yaml)
library(gitcreds)

### Functions modified from Baldan et al 2018 R package for aquaculture. 
### https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195732
### https://github.com/cran/RAC/tree/master/R

# Load required packages
library(qs)
library(qs2) 
library(terra)
library(readxl)
library(matrixStats)
library(furrr)
library(future)
library(dplyr)

make_label <- function(lab){lab %>% str_remove_all("_stat") %>% str_replace_all("_", " ") %>% str_to_title()}
fixnum <- function(n, digits = 4) {
  vapply(n, function(x) {
    str_flatten(c(rep("0", digits-nchar(as.character(x))), as.character(x)))
  }, character(1))
}
meanna <- function(x, ...) mean(x, na.rm = TRUE, ...)
minna <- function(x, ...) min(x, na.rm = TRUE, ...)
maxna <- function(x, ...) max(x, na.rm = TRUE, ...)
sdna <- function(x, ...) sd(x, na.rm = TRUE, ...)
sumna <- function(x, ...) sum(x, na.rm = TRUE, ...)
medianna <- function(x, ...) median(x, na.rm = TRUE, ...)
find_read <- function(path, pattern){
  file <- list.files(path, full.names = T, pattern = pattern)
  if (length(file) > 1) {print("Multiple files found - try again")} else {
    if (str_detect(file, ".qs")) {return(qs::qread(file))}
    if (str_detect(file, ".parquet")) {return(arrow::read_parquet(file))}
  }
}

# Parameters definition
# species_params['alpha']         [-] Feeding catabolism coefficient
# species_params['betaprot']      [-] Assimilation coefficient for protein - SUPERCEEDED by digestibility coefficient
# species_params['betalip']       [-] Assimilation coefficient for lipid - SUPERCEEDED by digestibility coefficient
# species_params['betacarb']      [-] Assimilation coefficient for carbohydrates - SUPERCEEDED by digestibility coefficient
# species_params['epsprot']       [J/gprot] Energy content of protein
# species_params['epslip']        [J/glip] Energy content of lipid
# species_params['epscarb']       [J/gcarb] Energy content of carbohydrate
# species_params['epsO2']         [J/gO2] Energy consumed by the respiration of 1g of oxygen
# species_params['pk']            [1/day] Temperature coefficient for the fasting catabolism
# species_params['k0']            [1/Celsius degree] Fasting catabolism at 0 Celsius degree
# species_params['m']             [-] Weight exponent for the anabolism
# species_params['n']             [-] Weight exponent for the catabolism
# species_params['betac']         [-] Shape coefficient for the H(Tw) function
# species_params['Tma']           [Celsius degree] Maximum lethal temperature  
# species_params['Toa']           [Celsius degree] Optimal temperature
# species_params['Taa']           [Celsius degree] Lowest feeding temperature
# species_params['omega']         [gO2/g] Oxygen consumption - weight loss ratio
# species_params['a']             [J/gtissue] Energy content of fish tissue
# species_params['k']             [-] Weight exponent for energy content
# species_params['eff']           [-] Food ingestion efficiency
# species_params['fcr']           [-] Food conversion ratio

get_farms <- function(farms_file, farm_ID, this_species){
  qread(farms_file) %>% 
    filter(model_name == this_species) %>% 
    select(-row_num) %>% 
    mutate(farm_id = row_number()) %>% 
    filter(farm_id == farm_ID)
}

get_feed_params <- function(file){
  df <- read.csv(file, header = F)
  values <- as.numeric(df$V1)
  names(values) <- df$V2
  values <- values[!is.na(values)]
  return(values)
}

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

feeding_rate <- function(water_temp, species_params) {
  exp(species_params['betac'] * (water_temp - species_params['Toa'])) * 
    ((species_params['Tma'] - water_temp)/(species_params['Tma'] - species_params['Toa']))^
    (species_params['betac'] * (species_params['Tma'] - species_params['Toa']))
}

food_prov_rate <- function(pop_params, water_temp, ing_pot, ing_pot_10, species_params) {
  # Use ifelse vectorization instead of individual if statements
  ifelse(
    water_temp > species_params['Taa'],
    ing_pot * (1 + rnorm(1, pop_params['overFmean'], pop_params['overFdelta'])),
    ing_pot_10
  ) # old formula: 0.25 * 0.066 * weight^0.75
}

# Apportion ingested feed into relevant components
app_feed <- function(provided, ingested, prop, macro, digestibility) {
  # Pre-compute common values and use vectorized operations
  provided_amount <- provided * prop * macro
  ingested_amount <- ingested * prop * macro
  assimilated <- ingested_amount * digestibility
  
  # Return only necessary values in a numeric vector
  c(provided = sum(provided_amount),
    ingested = sum(ingested_amount),
    uneaten = sum(provided_amount - ingested_amount),
    assimilated = sum(assimilated),
    excreted = sum(ingested_amount - assimilated))
}

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
    # Temperature response and feeding calculations
    result[i, 'rel_feeding'] <- feeding_rate(result[i, 'water_temp'], species_params)
    result[i, 'ing_pot'] <- ingmax * (result[i, 'weight']^species_params['m']) * result[i, 'rel_feeding']
    
    # Food provision and ingestion
    result[i, 'food_prov'] <- food_prov_rate(
      pop_params = pop_params, 
      water_temp = result[i, 'water_temp'],
      ing_pot = result[i, 'ing_pot'],
      ing_pot_10 = ingmax * (result[i, 'weight']^species_params['m']) * 0.1,
      species_params
    )
    result[i, 'food_enc'] <- species_params['eff'] * result[i, 'food_prov']
    result[i, 'ing_act'] <- min(result[i, 'food_enc'], result[i, 'ing_pot'])
    
    # Energy calculations
    result[i, 'E_somat'] <- species_params['a'] * result[i, 'weight']^species_params['k']
    
    # Process feed components - vectorized operations
    app_carbs <- app_feed(result[i, 'food_prov'], result[i, 'ing_act'],
                          feed_params[['Carbohydrates']]$proportion,
                          feed_params[['Carbohydrates']]$macro,
                          feed_params[['Carbohydrates']]$digest)
    app_lipids <- app_feed(result[i, 'food_prov'], result[i, 'ing_act'],
                           feed_params[['Lipids']]$proportion,
                           feed_params[['Lipids']]$macro,
                           feed_params[['Lipids']]$digest)
    app_proteins <- app_feed(result[i, 'food_prov'], result[i, 'ing_act'],
                             feed_params[['Proteins']]$proportion,
                             feed_params[['Proteins']]$macro,
                             feed_params[['Proteins']]$digest)
    
    # Store excretion and waste values
    result[i, c('C_excr', 'L_excr', 'P_excr')] <- c(app_carbs['excreted'], 
                                                    app_lipids['excreted'], 
                                                    app_proteins['excreted'])
    result[i, c('C_uneat', 'L_uneat', 'P_uneat')] <- c(app_carbs['uneaten'], 
                                                       app_lipids['uneaten'], 
                                                       app_proteins['uneaten'])
    
    # Energy assimilation
    result[i, 'E_assim'] <- app_carbs['assimilated'] * species_params['epscarb'] +
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

farm_growth <- function(pop_params, species_params, feed_params, water_temp, times, N_pop, nruns){
    
  days <- (times['t_start']:times['t_end'])*times['dt']
  
  # Generate all random values upfront
  init_weights <- rnorm(nruns, mean = pop_params['meanW'], sd = pop_params['deltaW'])
  ingmaxes <- rnorm(nruns, mean = pop_params['meanImax'], sd = pop_params['deltaImax'])
  
  # Run parallel simulation for individuals
  mc_results <- furrr::future_map2(init_weights, ingmaxes, function(init_w, ing_m) {
    mat <- fish_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = water_temp,
      feed_params = feed_params,
      times = times,
      init_weight = init_w,
      ingmax = ing_m
    ) %>% unname()
  }, .options = furrr::furrr_options(seed = TRUE))
  
  stat_names <- c("days", "weight", "dw", "water_temp", "T_response", "P_excr", "L_excr", "C_excr", "P_uneat", 
                  "L_uneat", "C_uneat", "food_prov", "food_enc", "rel_feeding", "ing_pot", "ing_act", "E_assim", 
                  "E_somat", "anab", "catab", "O2", "NH4")

  # Consolidate all individuals into a farm (with population = nruns)
  all_results <- lapply(1:length(stat_names), function(col_idx) {
    t(sapply(mc_results, function(mat_idx) {
      mat_idx[, col_idx]
    }))
  }) %>% setNames(stat_names)
  
  # Some stats need to be summed/added
  all_results[["total_excr"]] <- all_results[["P_excr"]] + all_results[["L_excr"]] + all_results[["C_excr"]]
  all_results[["total_uneat"]] <- all_results[["P_uneat"]] + all_results[["L_uneat"]] + all_results[["C_uneat"]]
  all_results[["metab"]] <- all_results[["anab"]] - all_results[["catab"]]
  all_results[["biomass"]] <- all_results[["weight"]]
  
  all_results <- lapply(2:length(all_results), function(col_idx) {
      cbind(colMeans(all_results[[col_idx]]), matrixStats::colSds(all_results[[col_idx]])) %>% 
      as.matrix() %>% unname()
  }) %>% setNames(names(all_results)[2:length(names(all_results))])

  # Some stats should be multiplied by the farm population (Npop)
  pop_names <- c("biomass", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "ing_act", 
                 "total_excr", "total_uneat", "O2", "NH4", "food_prov")
  for (stat_nm in pop_names) {
    all_results[[stat_nm]][,1] <- all_results[[stat_nm]][,1] * N_pop[1:length(days)]
    all_results[[stat_nm]][,2] <- all_results[[stat_nm]][,2] * N_pop[1:length(days)]
  }
  
  out_list <- lapply(1:length(all_results), function(col_idx) {
    cbind(days, all_results[[col_idx]]) %>% 
      as.matrix() %>% unname() 
  }) %>% setNames(paste0(names(all_results), "_stat"))
  
  return(out_list)
}

# This is identical to the farm_growth function except without the Monte-Carlo sampling of initial weights (all uniform)
uni_farm_growth <- function(pop_params, species_params, feed_params, water_temp, times, N_pop){
  
  days <- (times['t_start']:times['t_end'])*times['dt']
  
  # Run parallel simulation for individuals
  mc_results2 <- fish_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = water_temp,
      feed_params = feed_params,
      times = times,
      init_weight = pop_params['meanW'],
      ingmax = pop_params['meanImax']
    ) %>% unname()
  
  stat_names <- c("days", "weight", "dw", "water_temp", "T_response", "P_excr", "L_excr", "C_excr", "P_uneat", 
                  "L_uneat", "C_uneat", "food_prov", "food_enc", "rel_feeding", "ing_pot", "ing_act", "E_assim", 
                  "E_somat", "anab", "catab", "O2", "NH4")
  
  all_results2 <- setNames(split(t(mc_results2), row(t(mc_results2))), stat_names)
  all_results2 <- all_results2[-1]

  # Some stats need to be summed/added
  all_results2[["total_excr"]] <- all_results2[["P_excr"]] + all_results2[["L_excr"]] + all_results2[["C_excr"]]
  all_results2[["total_uneat"]] <- all_results2[["P_uneat"]] + all_results2[["L_uneat"]] + all_results2[["C_uneat"]]
  all_results2[["metab"]] <- all_results2[["anab"]] - all_results2[["catab"]]
  all_results2[["biomass"]] <- all_results2[["weight"]]
  
  # Some stats should be multiplied by the farm population (Npop)
  pop_names <- c("biomass", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "ing_act", 
                 "total_excr", "total_uneat", "O2", "NH4", "food_prov")
  for (stat_nm in pop_names) {
    all_results2[[stat_nm]] <- all_results2[[stat_nm]] * N_pop[1:length(days)]
  }
  
  out_list <- lapply(1:length(all_results2), function(col_idx) {
    cbind(days, all_results2[[col_idx]]) %>% 
      as.matrix() %>% unname() 
  }) %>% setNames(paste0(names(all_results2), "_stat"))
  
  return(out_list)
}

cite_pack <- function(package_name) {
  version <- packageVersion(package_name)
  paste0("package `", package_name, "`, version ", version, " [@", package_name, "]")
}

cite_packages <- function(packages) {
  # Generate citations for each package
  citations <- sapply(packages, function(pkg) {
    version <- packageVersion(pkg)
    paste0("`", pkg, "`, version ", version, " [@", pkg, "]")
  })
  
  # Handle different cases based on number of packages
  if (length(packages) == 1) {
    return(paste0("package ", citations[1]))
  } else if (length(packages) == 2) {
    return(paste0("packages ", citations[1], " and ", citations[2]))
  } else {
    # More than 2 packages: use commas and "and" before last
    first_part <- paste(citations[-length(citations)], collapse = ", ")
    last_part <- citations[length(citations)]
    return(paste0("packages ", first_part, ", and ", last_part))
  }
}

create_boxes <- function(box_list, label_positions, offset_deg) {
  boxes <- list()
  labels <- list()
  
  for(i in seq_along(box_list)) {
    name <- names(box_list)[i]
    coords <- box_list[[i]]
    pos <- label_positions[i]
    
    # Create polygon from bounding box coordinates
    box_coords <- matrix(c(
      coords[1], coords[3],  # lonmin, latmin
      coords[2], coords[3],  # lonmax, latmin  
      coords[2], coords[4],  # lonmax, latmax
      coords[1], coords[4],  # lonmin, latmax
      coords[1], coords[3]   # close polygon
    ), ncol = 2, byrow = TRUE)
    
    # Create sf polygon
    poly <- st_polygon(list(box_coords))
    boxes[[i]] <- st_sfc(poly, crs = 4326)
    
    # Calculate label position based on specification
    label_coords <- get_label_position(coords, pos, offset_deg[i])
    label_point <- st_point(c(label_coords[1], label_coords[2]))
    labels[[i]] <- st_sfc(label_point, crs = 4326)
  }
  
  # Combine into sf objects
  boxes_sf <- st_sf(
    name = names(box_list),
    letter = LETTERS[1:length(box_list)],
    geometry = do.call(c, boxes)
  )
  
  labels_sf <- st_sf(
    name = names(box_list), 
    letter = LETTERS[1:length(box_list)],
    geometry = do.call(c, labels)
  )
  
  return(list(boxes = boxes_sf, labels = labels_sf))
}

# Helper function to calculate label positions
get_label_position <- function(coords, position, offset) {
  lonmin <- coords[1]
  lonmax <- coords[2]
  latmin <- coords[3]
  latmax <- coords[4]
  
  # Parse position string
  parts <- strsplit(position, "_")[[1]]
  vertical <- parts[1]    # "top" or "bottom"
  horizontal <- parts[2]  # "left" or "right"
  location <- parts[3]    # "inside" or "outside"
  
  # Base coordinates for corners
  if(vertical == "bottom") {base_lat <- latmin} else {base_lat <- latmax}
  if(horizontal == "right") {base_lon <- lonmax} else {base_lon <- lonmin}
  
  # Apply offset for outside positioning
  if(location == "outside") {
    if(horizontal == "left") {base_lon <- base_lon-offset} else {base_lon <- base_lon+offset}
    if(vertical == "top") {base_lat <- base_lat + offset} else {base_lat <- base_lat - offset}
  }
  
  return(c(base_lon, base_lat))
}

lim_robin <- function(lons = c(144.0, 149.5), lats = c(-39.75, -44.00)) {
  coords <- data.frame(lon = lons, lat = lats) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(crs = "+proj=robin") %>% 
    st_coordinates()
  list(
    xlims = range(coords[, "X"]),
    ylims = range(coords[, "Y"])
  )
}

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

# nolint end
