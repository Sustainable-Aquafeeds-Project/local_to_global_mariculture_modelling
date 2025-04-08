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

fixnum <- function(n, digits = 4) {
  str_flatten(c(rep("0", digits-nchar(as.character(n))), as.character(n)))
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

feeding_rate <- function(water_temp, species_params){
  exp(species_params['betac']*(water_temp-species_params['Toa'])) * 
    ((species_params['Tma']-water_temp)/(species_params['Tma']-species_params['Toa']))^(species_params['betac']*(species_params['Tma']-species_params['Toa']))
}

food_prov_rate <- function(rel_feeding, ing_pot, ing_pot_10){
  ifelse(
    rel_feeding > 0.1, 
    yes = ing_pot + rnorm(1, 0.035, 0.0075)*ing_pot, 
    no = ing_pot_10 # no = 0.25 * 0.066 * weight^0.75
  )
}

app_feed <- function(provided, ingested, prop, macro, digestibility) {
  app <- matrix(0, nrow = length(prop), ncol = 5, dimnames = list(NULL, c('provided', 'ingested', 'uneaten', 'assimilated', 'excreted')))
  app[,'provided'] <- provided * prop * macro
  app[,'ingested'] <- ingested * prop * macro
  app[,'assimilated'] <- app[,'ingested'] * digestibility
  app[,'uneaten'] <- app[,'provided'] - app[,'ingested']
  app[,'excreted'] <- app[,'ingested'] - app[,'assimilated']
  return(colSums(app))
}

fish_growth <- function(pop_params, species_params, water_temp, feed_params, times, init_weight, ingmax){
  # Weight and resource
  days <- (times['t_start']:times['t_end'])*times['dt']
  
  # vectors of the other important growth metrics: intake, functional response, assimilated energy and growth increment
  rel_feeding <- T_response <- ing_pot <- food_prov <- food_enc <- ing_act <- E_somat <- E_assim <- P_excr <- L_excr <- C_excr <- P_uneat <- L_uneat <- C_uneat <- anab <- catab <- O2 <- NH4 <- dw <- weight <- rep(0, length(days))

  # start size in grams
  weight[1] <- init_weight
  
  # EQUATIONS
  for (i in 1:(length(days-1))){
    # Temperature dependence for relative maximum ingestion rate
    rel_feeding[i] <- feeding_rate(water_temp[i], species_params)
    
    # Ingested mass [g/d]
    ing_pot[i] <- ingmax * (weight[i]^species_params['m']) * rel_feeding[i]
    food_prov[i] <- food_prov_rate(rel_feeding[i], ing_pot[i], ing_pot_10 = ingmax * (weight[i]^species_params['m']) * 0.1)
    food_enc[i] <- species_params['eff'] * food_prov[i] # encounter rate
    ing_act[i] <- min(food_enc[i], ing_pot[i]) # [g/d] Actual ingestion rate cannot be more than encountered food
    
    # Energy content of somatic tissue [J/g] Source: Lupatsch et al. (2003)
    E_somat[i] = species_params['a'] * weight[i]^species_params['k']
    
    # Apportion food to ingested/uneaten, assimilated/excreted [g/d]
    app_carbs <- app_feed(food_prov[i], ing_act[i], feed_params[['Carbohydrates']]$proportion, feed_params[['Carbohydrates']]$macro, feed_params[['Carbohydrates']]$digest)
    app_lipids <- app_feed(food_prov[i], ing_act[i], feed_params[['Lipids']]$proportion, feed_params[['Lipids']]$macro, feed_params[['Lipids']]$digest)
    app_proteins <- app_feed(food_prov[i], ing_act[i], feed_params[['Proteins']]$proportion, feed_params[['Proteins']]$macro, feed_params[['Proteins']]$digest)
    
    # [J/d] Assimilated energy
    E_assim[i] <- app_carbs['assimilated']*species_params['epscarb'] + app_lipids['assimilated']*species_params['epslip'] + app_proteins['assimilated']*species_params['epsprot']
    
    # Excretion (faeces)
    C_excr[i] <- app_carbs['excreted']     # Excreted carbohydrates [g/d]
    L_excr[i] <- app_lipids['excreted']    # Excreted lipids [g/d]
    P_excr[i] <- app_proteins['excreted']  # Excreted proteins [g/d]
    
    # Uneaten feed
    C_uneat[i] <- app_carbs['uneaten']       # Carbohydrates to uneaten feed [g/d]
    L_uneat[i] <- app_lipids['uneaten']      # Lipids to uneaten feed [g/d]
    P_uneat[i] <- app_proteins['uneaten']    # Proteins to uneaten feed [g/d]
    
    # Exponential Temperature dependence for catabolism
    T_response[i] <- exp(species_params['pk']*water_temp[i])
    
    # Metabolism terms
    anab[i]  <- E_assim[i]*(1-species_params['alpha'])                    # Net anabolism [J/d]
    catab[i] <- species_params['epsO2'] * species_params['k0'] * T_response[i] * (weight[i]^species_params['n'])*species_params['omega']    # Fasting catabolism [J/d]
    
    # O2 and NH4 produced
    O2[i]  <- catab[i]/species_params['epsO2']        # O2 consumed [g02/d]
    NH4[i] <- O2[i]*0.06                              # NH4 produced [gN/d]
    
    # Mass balance
    dw[i]  <- (anab[i] - catab[i]) / E_somat[i]       # weight increment [g/d]
    
    weight[i+1] <- weight[i] + dw[i]*times['dt']
  }
  
  # Function outputs
  output <- cbind(days, weight = weight[1:length(days)], dw, water_temp = water_temp[1:length(days)], T_response, P_excr, L_excr, C_excr, P_uneat, L_uneat, C_uneat, food_prov, food_enc, rel_feeding, ing_pot, ing_act, E_assim, E_somat, anab, catab, O2, NH4)
  
  return(output)
}

farm_growth <- function(pop_params, species_params, feed_params, water_temp, times, N_pop, nruns){

  days <- (times['t_start']:times['t_end'])*times['dt']
  
  # Initiate matrices to fill for each population iteration
  weight_mat <- biomass_mat <- dw_mat <- SGR_mat <- E_somat_mat <- P_excr_mat <-  L_excr_mat <- C_excr_mat <- P_uneat_mat <- L_uneat_mat <- C_uneat_mat <- ing_act_mat <- anab_mat <- catab_mat <- O2_mat <- NH4_mat <- food_prov_mat <- rel_feeding_mat <- T_response_mat <- total_excr_mat <- total_uneat_mat <- matrix(data = 0, nrow = nruns, ncol = length(days)) 
  
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
    # if(n %in% seq(500,5000,500)){print(paste(n, "runs done at", Sys.time()))}
    
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
    O2_stat = cbind(colMeans(O2_mat), colSds(O2_mat)),
    NH4_stat = cbind(colMeans(NH4_mat), colSds(NH4_mat)),
    food_prov_stat = cbind(colMeans(food_prov_mat), colSds(food_prov_mat)),
    rel_feeding_stat = cbind(colMeans(rel_feeding_mat), colSds(rel_feeding_mat)),
    T_response_stat = cbind(colMeans(T_response_mat), colSds(T_response_mat)),
    total_excr_mat = cbind(colMeans(total_excr_mat), colSds(total_excr_mat)),
    total_uneat_mat = cbind(colMeans(total_uneat_mat), colSds(total_uneat_mat))
  )
  return(out_list)
}

days_to_CS <- function(weight_stat, days, CS){
  days[weight_stat[which(weight_stat > CS)]][1]
}

post_process <- function(path, farm_ID, Feedfype, out_loop, times, N_pop, CS) {
  times['t_start']=times[1]         # Integration beginning
  times['t_end']=times[2]           # Integration end
  
  # Extracts outputs from the output list
  weight_stat=out_loop[[1]]
  Bio_stat = out_loop[[2]]
  Dw_stat = out_loop[[3]]
  E_somat_stat = out_loop[[4]]
  P_excr_stat=out_loop[[5]]
  L_excr_stat=out_loop[[6]]
  C_excr_stat=out_loop[[7]]
  P_uneat_stat=out_loop[[8]]
  L_uneat_stat=out_loop[[9]]
  C_uneat_stat=out_loop[[10]]
  ing_act_stat=out_loop[[11]]
  A_stat=out_loop[[12]]
  C_stat=out_loop[[13]]
  NH4_stat=out_loop[[14]]
  O2_stat=out_loop[[15]]
  Resource_stat = out_loop[[16]]
  rel_feeding_stat=out_loop[[17]]
  T_response_stat=out_loop[[18]]
  
  
  # Days to commercial size
  # Lower bound
  foo <- function(w, S) {
    which(w > S)[1]
  }
  arg <- as.data.frame(weight_stat[, 1] - weight_stat[, 2]) # Mean - SD
  days <- apply(arg, 1, foo, S = CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex) == 0) {
    Lb_daysToSize = "Not reaching the commercial size"
  } else{
    Lb_daysToSize <- min(NonNAindex)
  }
  
  # Mean
  foo <- function(w, S) {
    which(w > S)[1]
  }
  arg = as.data.frame(weight_stat[, 1])
  days <- apply(arg, 1, foo, S = CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex) == 0) {
    Mean_daysToSize = "Not reaching the commercial size"
  } else{
    Mean_daysToSize <- min(NonNAindex)
  }
  
  # Upper bound
  foo <- function(w, S) {
    which(w > S)[1]
  }
  arg = as.data.frame(weight_stat[, 1] + weight_stat[, 2])
  days <- apply(arg, 1, foo, S = CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex) == 0) {
    Ub_daysToSize = "Not reaching the commercial size"
  } else{
    Ub_daysToSize <- min(NonNAindex)
  }
  
  # List containing days to size
  daysToSize<-as.list(cbind(Ub_daysToSize,Mean_daysToSize,Lb_daysToSize))
  
  
  
  out_post = list(weight_stat, Bio_stat, Dw_stat, E_somat_stat, P_excr_stat, L_excr_stat, C_excr_stat, P_uneat_stat, L_uneat_stat, C_uneat_stat, A_stat, C_stat, NH4_stat,O2_stat,  Resource_stat, rel_feeding_stat, T_response_stat, daysToSize, N_pop)
  
 
  
  
  # PLOT RESULTS
  #shorter and longer version of days due to lag being used in some functions
  
  days <- seq(from = times['t_start'], to = times['t_end']-1, by = 1) # create a dates vector to plot results
  days2 <- seq(times['t_start'], by = 1, length = times['t_end']-times['t_start']+1) # create a dates vector to plot results
  
  
  
  
  
  # Plot weight
  
  weight_df = data.frame(days = days2, weight = weight_stat[,1], lower_bound = weight_stat[,1]-weight_stat[,2], upper_bound = weight_stat[,1]+weight_stat[,2])
  
  #plot weight    
  ggplot(data = weight_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = weight))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Weight (g)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/weight_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  
  # plot biomass
  
  
  biomass_df = data.frame(days = days2, biomass = Bio_stat[,1], lower_bound = Bio_stat[,1]-Bio_stat[,2], upper_bound = Bio_stat[,1]+Bio_stat[,2])
  
  
  ggplot(data = biomass_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = biomass))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Total farm biomass (tonnes)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = function(label){label/1e+6})
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/biomass_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  # plot dW
  
  
  dw_df = data.frame(days = days, dw = Dw_stat[,1][-730], lower_bound = Dw_stat[,1][-730]-Dw_stat[,2][-730], upper_bound = Dw_stat[,1][-730]+Dw_stat[,2][-730])
  
  
  ggplot(data = dw_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = dw))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = bquote(Delta~ Body~weight~g~day^1))+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/dw_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  
  #plot energy content of somatic tissue
  
  E_somat_df = data.frame(days = days, E_somat = E_somat_stat[,1][-730], lower_bound = E_somat_stat[,1][-730]-E_somat_stat[,2][-730], upper_bound = E_somat_stat[,1][-730]+E_somat_stat[,2][-730])
  
  
  ggplot(data = E_somat_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = E_somat))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Energy content of somatic tissue (J/g)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))
  
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/E_somat_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  # plot excretion
  
  excretion_df <- 
    rbind(
      data.frame(days = days2, 
                 mean_excretion = P_excr_stat[,1], 
                 lower_bound = P_excr_stat[,1]-P_excr_stat[,2], 
                 upper_bound = P_excr_stat[,1]+P_excr_stat[,2],
                 nutrient = "Protein"),
      data.frame(days = days2, 
                 mean_excretion = L_excr_stat[,1], 
                 lower_bound = L_excr_stat[,1]-L_excr_stat[,2], 
                 upper_bound = L_excr_stat[,1]+L_excr_stat[,2],
                 nutrient = "Lipid"),
      data.frame(days = days2, 
                 mean_excretion = C_excr_stat[,1], 
                 lower_bound = C_excr_stat[,1]-C_excr_stat[,2], 
                 upper_bound = C_excr_stat[,1]+C_excr_stat[,2],
                 nutrient = "Carbohydrates")
    )
  
  
  ggplot(data = excretion_df |> filter(days != 730))+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound, fill = nutrient), alpha = 0.2)+
    geom_line(aes(x = days, y = mean_excretion, colour = nutrient))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Excretion (kg/day)")+
    theme(text=element_text(size=8),
          legend.position = "inside",
          legend.title = element_blank(),
          legend.position.inside= c(0.2, 0.8),
          legend.background = element_rect(fill="transparent"))+
    scale_y_continuous(labels = function(label){label/1e+3})
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/excretion_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  # plot wasted feed
  
  feed_waste_df <- 
    rbind(
      data.frame(days = days2, 
                 mean_waste = P_uneat_stat[,1], 
                 lower_bound = P_uneat_stat[,1]-P_uneat_stat[,2], 
                 upper_bound = P_uneat_stat[,1]+P_uneat_stat[,2],
                 nutrient = "Protein"),
      data.frame(days = days2, 
                 mean_waste = L_uneat_stat[,1], 
                 lower_bound = L_uneat_stat[,1]-L_uneat_stat[,2], 
                 upper_bound = L_uneat_stat[,1]+L_uneat_stat[,2],
                 nutrient = "Lipid"),
      data.frame(days = days2, 
                 mean_waste = C_uneat_stat[,1], 
                 lower_bound = C_uneat_stat[,1]-C_uneat_stat[,2], 
                 upper_bound = C_uneat_stat[,1]+C_uneat_stat[,2],
                 nutrient = "Carbohydrates")
    )
  
  ggplot(data = feed_waste_df |> filter(days != 730))+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound, fill = nutrient), alpha = 0.2)+
    geom_line(aes(x = days, y = mean_waste, colour = nutrient))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Feed waste (kg/day)")+
    theme(text=element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.25, 0.8),
          legend.background = element_rect(fill="transparent"),
          legend.title = element_blank())+
    scale_y_continuous(labels = function(label){label/1e+3})
  
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/feed_waste_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  #plot ingestion
  
  
  ing_act_df = data.frame(days = days, ing_act = ing_act_stat[,1][-730], lower_bound = ing_act_stat[,1][-730]-ing_act_stat[,2][-730], upper_bound = ing_act_stat[,1][-730]+ing_act_stat[,2][-730])
  
  
  
  ggplot(data = ing_act_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = ing_act))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Ingestion (kg/day)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = function(label){label/1e+3})
  
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/actual_ingestion_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  metab_df =  rbind(
    data.frame(days = days, value = A_stat[,1][-730], lower_bound = A_stat[,1][-730]-A_stat[,2][-730], upper_bound = A_stat[,1][-730]+A_stat[,2][-730], metabolic_f = "Anabolism"),
    data.frame(days = days, value = C_stat[,1][-730], lower_bound = C_stat[,1][-730]-C_stat[,2][-730], upper_bound = C_stat[,1][-730]+C_stat[,2][-730], metabolic_f = "Catabolism")
  )
  
  ggplot(data = metab_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound, fill = metabolic_f), alpha = 0.2)+
    geom_line(aes(x = days, y = value, colour = metabolic_f))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Net matbolic rate (J/day)")+
    theme(text=element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.25, 0.8),
          legend.background = element_rect(fill="transparent"),
          legend.title = element_blank())+
    scale_y_continuous(labels = function(label){label/1e+3})
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/metabolic_rate_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  #Plot O2 consumption
  
  O2_df = data.frame(days = days, o2_consumption = O2_stat[,1][-730], lower_bound = O2_stat[,1][-730]-O2_stat[,2][-730], upper_bound = O2_stat[,1][-730]+O2_stat[,2][-730])
  
  ggplot(data = O2_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = o2_consumption))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Oxygen consumption (kg/day)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = function(labels) labels/1e+3)
  
  
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/o2_consumption_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  #plot ammonia
  
  NH4_df = data.frame(days = days, nh4_production = NH4_stat[,1][-730], lower_bound = NH4_stat[,1][-730]-NH4_stat[,2][-730], upper_bound = NH4_stat[,1][-730]+NH4_stat[,2][-730])
  
  
  ggplot(data = NH4_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = nh4_production))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Ammonia production (kg N/ day)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = function(labels) labels/1e+3)
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/nh4_production_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  food_prov_df = data.frame(days = days, feed_resource = Resource_stat[,1][-730], lower_bound = Resource_stat[,1][-730]-Resource_stat[,2][-730], upper_bound = Resource_stat[,1][-730]+Resource_stat[,2][-730])
  
  
  ggplot(data = resource_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = feed_resource))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Available feed (kg/day)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = function(labels) labels/1e+3, limits = c(0,9000000))
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/available_feed_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  
  
  #plot temperature functions
  
  temp_f_df =  rbind(
    data.frame(days = days, value = rel_feeding_stat[,1][-730], lower_bound = rel_feeding_stat[,1][-730]-rel_feeding_stat[,2][-730], upper_bound = rel_feeding_stat[,1][-730]+rel_feeding_stat[,2][-730], metabolic_f = "Anabolism"),
    data.frame(days = days, value = T_response_stat[,1][-730], lower_bound = T_response_stat[,1][-730]-T_response_stat[,2][-730], upper_bound = T_response_stat[,1][-730]+T_response_stat[,2][-730], metabolic_f = "Catabolism")
  )
  
  
  
  ggplot(data = temp_f_df)+
    #geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound, fill = metabolic_f), alpha = 0.2)+
    geom_line(aes(x = days, y = value, colour = metabolic_f))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Temperature dependence function")+
    theme(text=element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.25, 0.5),
          legend.background = element_rect(fill="transparent"),
          legend.title = element_blank())+
    scale_y_continuous(labels = function(label){label/1e+3})
  
  ggsave(filename = file.path(path, sprintf("figures/outputs/%s/temp_function_%s.jpeg", Feed_type, farm_ID)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  #save data
  
  qsave(x= weight_df, file = file.path(path, sprintf("data_products/model_outputs/%s/weight_output_%s.qs", Feed_type, farm_ID)))
  qsave(x = biomass_df, file = file.path(path, sprintf("data_products/model_outputs/%s/biomass_output_%s.qs", Feed_type, farm_ID)))
  qsave(x = dw_df, file = file.path(path, sprintf("data_products/model_outputs/%s/dw_output_%s.qs", Feed_type, farm_ID)))
  qsave(E_somat_df, file = file.path(path, sprintf("data_products/model_outputs/%s/E_somat_output_%s.qs", Feed_type, farm_ID)))
  qsave(excretion_df, file = file.path(path, sprintf("data_products/model_outputs/%s/excretion_output_%s.qs", Feed_type, farm_ID)))
  qsave(feed_waste_df, file = file.path(path, sprintf("data_products/model_outputs/%s/feed_waste_output_%s.qs", Feed_type, farm_ID)))
  qsave(ing_act_df, file = file.path(path, sprintf("data_products/model_outputs/%s/ing_act_output_%s.qs", Feed_type, farm_ID)))
  qsave(metab_df, file = file.path(path, sprintf("data_products/model_outputs/%s/metabolic_rate_output_%s.qs", Feed_type, farm_ID)))
  qsave(O2_df, file = file.path(path, sprintf("data_products/model_outputs/%s/O2_consumption_output_%s.qs", Feed_type, farm_ID)))
  qsave(NH4_df, file = file.path(path, sprintf("data_products/model_outputs/%s/NH4_production_output_%s.qs", Feed_type, farm_ID)))
  qsave(resource_df, file = file.path(path, sprintf("data_products/model_outputs/%s/feed_available_output_%s.qs", Feed_type, farm_ID)))
  qsave(temp_f_df, file = file.path(path, sprintf("data_products/model_outputs/%s/temp_function_output_%s.qs", Feed_type, farm_ID)))
  qsave(daysToSize, file = file.path(path, sprintf("data_products/model_outputs/%s/days_to_size_output_%s.qs", Feed_type, farm_ID)))

}

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

# Function that solves the population dynamics equations including discontinuities
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


