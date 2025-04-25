### Functions modified from Baldan et al 2018 R package for aquaculture. 
### https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195732
### https://github.com/cran/RAC/tree/master/R

# Load required packages silently
suppressPackageStartupMessages({
  library(qs)
  library(qs2) 
  library(terra)
  library(readxl)
})

fixnum <- function(n, digits = 4) {str_flatten(c(rep("0", digits-nchar(as.character(n))), as.character(n)))}
meanna <- function(x, ...) mean(x, na.rm = TRUE, ...)
minna <- function(x, ...) min(x, na.rm = TRUE, ...)
maxna <- function(x, ...) max(x, na.rm = TRUE, ...)
sdna <- function(x, ...) sd(x, na.rm = TRUE, ...)
sumna <- function(x, ...) sum(x, na.rm = TRUE, ...)
medianna <- function(x, ...) median(x, na.rm = TRUE, ...)

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

food_prov_rate <- function(pop_params, rel_feeding, ing_pot, ing_pot_10) {
    # Use ifelse vectorization instead of individual if statements
    ifelse(
      rel_feeding > 0.1,
      ing_pot * (1 + rnorm(1, pop_params['overFmean'], pop_params['overFdelta'])),
      ing_pot_10
      ) # old formula: 0.25 * 0.066 * weight^0.75
}

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
          rel_feeding = result[i, 'rel_feeding'],
          ing_pot = result[i, 'ing_pot'],
          ing_pot_10 = ingmax * (result[i, 'weight']^species_params['m']) * 0.1
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
  
  # Initiate matrices to fill for each population iteration
  weight_mat <- biomass_mat <- dw_mat <- SGR_mat <- E_somat_mat <- P_excr_mat <-  L_excr_mat <- C_excr_mat <- P_uneat_mat <- L_uneat_mat <- C_uneat_mat <- ing_act_mat <- anab_mat <- catab_mat <- O2_mat <- NH4_mat <- food_prov_mat <- rel_feeding_mat <- T_response_mat <- total_excr_mat <- total_uneat_mat <- 
    matrix(data = 0, nrow = nruns, ncol = length(days)) 
  
  init_weight <- rnorm(nruns, mean = pop_params['meanW'], sd = pop_params['deltaW'])
  ingmax <- rnorm(nruns, mean = pop_params['meanImax'], sd = pop_params['deltaImax'])
  
  for(n in 1:nruns){
    ind_output <- fish_growth(
      pop_params = pop_params,
      species_params = species_params,
      water_temp = water_temp,
      feed_params = feed_params,
      times = times,
      init_weight = init_weight[n],
      ingmax = ingmax[n]
    )
    # Append to matrix
    weight_mat[n,]      <- ind_output[,'weight']
    biomass_mat[n,]     <- ind_output[,'weight']*N_pop[1:length(days)]
    dw_mat[n,]          <- ind_output[,'dw']
    SGR_mat[n,]         <- 100 * (exp((log(weight_mat[n,])-log(weight_mat[n,1]))/(ind_output[,'days'])) - 1)
    E_somat_mat[n,]     <- ind_output[,'E_somat']
    P_excr_mat[n,]      <- ind_output[,'P_excr']*N_pop[1:length(days)]
    L_excr_mat[n,]      <- ind_output[,'L_excr']*N_pop[1:length(days)]
    C_excr_mat[n,]      <- ind_output[,'C_excr']*N_pop[1:length(days)]
    P_uneat_mat[n,]     <- ind_output[,'P_uneat']*N_pop[1:length(days)]
    L_uneat_mat[n,]     <- ind_output[,'L_uneat']*N_pop[1:length(days)]
    C_uneat_mat[n,]     <- ind_output[,'C_uneat']*N_pop[1:length(days)]
    ing_act_mat[n,]     <- ind_output[,'ing_act']*N_pop[1:length(days)]
    anab_mat[n,]        <- ind_output[,'anab']
    catab_mat[n,]       <- ind_output[,'catab']
    O2_mat[n,]          <- ind_output[,'O2']
    NH4_mat[n,]         <- ind_output[,'NH4']*N_pop[1:length(days)]
    food_prov_mat[n,]   <- ind_output[,'food_prov']*N_pop[1:length(days)]
    rel_feeding_mat[n,] <- ind_output[,'rel_feeding']
    T_response_mat[n,]  <- ind_output[,'T_response']
    total_excr_mat[n,]  <- (ind_output[,'P_excr'] + ind_output[,'L_excr'] + ind_output[,'C_excr']) * N_pop[1:length(days)]
    total_uneat_mat[n,] <- (ind_output[,'P_uneat'] + ind_output[,'L_uneat'] + ind_output[,'C_uneat']) * N_pop[1:length(days)]
  }

  out_list <- list(
      weight_stat = cbind(colMeans(weight_mat), colSds(weight_mat)),
      biomass_stat = cbind(colMeans(biomass_mat), colSds(biomass_mat)),
      dw_stat = cbind(colMeans(dw_mat), colSds(dw_mat)),
      SGR_stat = cbind(colMeans(SGR_mat), colSds(SGR_mat)),
      E_somat_stat = cbind(colMeans(E_somat_mat), colSds(E_somat_mat)),
      P_excr_stat = cbind(colMeans(P_excr_mat), colSds(P_excr_mat)),
      L_excr_stat = cbind(colMeans(L_excr_mat), colSds(L_excr_mat)),
      C_excr_stat = cbind(colMeans(C_excr_mat), colSds(C_excr_mat)),
      P_uneat_stat = cbind(colMeans(P_uneat_mat), colSds(P_uneat_mat)),
      L_uneat_stat = cbind(colMeans(L_uneat_mat), colSds(L_uneat_mat)),
      C_uneat_stat = cbind(colMeans(C_uneat_mat), colSds(C_uneat_mat)),
      ing_act_stat = cbind(colMeans(ing_act_mat), colSds(ing_act_mat)),
      anab_stat = cbind(colMeans(anab_mat), colSds(anab_mat)),
      catab_stat = cbind(colMeans(catab_mat), colSds(catab_mat)),
      NH4_stat = cbind(colMeans(NH4_mat), colSds(NH4_mat)),
      O2_stat = cbind(colMeans(O2_mat), colSds(O2_mat)),
      food_prov_stat = cbind(colMeans(food_prov_mat), colSds(food_prov_mat)),
      rel_feeding_stat = cbind(colMeans(rel_feeding_mat), colSds(rel_feeding_mat)),
      T_response_stat = cbind(colMeans(T_response_mat), colSds(T_response_mat))
  )
  return(out_list)
}

