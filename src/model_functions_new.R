### Main model functions here are modified from Baldan et al 2018 R package for aquaculture.
### https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195732
### https://github.com/cran/RAC/tree/master/R

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
library(here)

# species_params named vector — parameter definitions
# --- Bioenergetic parameters ---
# species_params['alpha']         [-]            Feeding catabolism coefficient
# species_params['epsprot']       [J/gprot]      Energy content of protein
# species_params['epslip']        [J/glip]       Energy content of lipid
# species_params['epscarb']       [J/gcarb]      Energy content of carbohydrate
# species_params['epsO2']         [J/gO2]        Energy consumed by the respiration of 1 g of oxygen
# species_params['pk']            [1/day]        Temperature coefficient for the fasting catabolism
# species_params['k0']            [1/°C]         Fasting catabolism at 0°C
# species_params['m']             [-]            Weight exponent for the anabolism
# species_params['n']             [-]            Weight exponent for the catabolism
# species_params['betac']         [-]            Shape coefficient for the H(Tw) function
# species_params['Tma']           [°C]           Maximum feeding temperature
# species_params['Toa']           [°C]           Optimal feeding temperature
# species_params['Taa']           [°C]           Lowest (minimum) feeding temperature
# species_params['omega']         [gO2/g]        Oxygen consumption–weight loss ratio
# species_params['a']             [J/gtissue]    Energy content of fish tissue
# species_params['k']             [-]            Weight exponent for energy content
# species_params['eff']           [-]            Food ingestion efficiency
# species_params['fcr']           [-]            Food conversion ratio
# --- Population / individual-variability parameters ---
# species_params['meanW']         [g]            Mean initial individual weight
# species_params['deltaW']        [g]            SD of initial individual weight
# species_params['meanImax']      [g/day]        Mean maximum ingestion rate
# species_params['deltaImax']     [g/day]        SD of maximum ingestion rate
# species_params['overFmean']     [-]            Mean overfeeding multiplier (added to 1)
# species_params['overFdelta']    [-]            SD of overfeeding multiplier


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

#' Convert protein mass to nitrogen mass
#'
#' Uses the standard 6.25 nitrogen-to-protein conversion factor (Jones factor).
#'
#' @param P Numeric. Protein mass (g).
#'
#' @return Numeric. Nitrogen mass (g).
get_nitrogen <- function(P) {unname(P / 6.25)}


#' Convert macronutrient masses to carbon mass
#'
#' Estimates the total carbon content (g) from the masses of lipid, carbohydrate
#' and protein, using fixed elemental conversion factors.
#'
#' @param P Numeric. Protein mass (g).
#' @param L Numeric. Lipid mass (g).
#' @param C Numeric. Carbohydrate mass (g).
#'
#' @return Numeric. Carbon mass (g).
get_carbon <- function(P, L, C) {
  carbon_1 <- L * 0.75
  carbon_2 <- C * 0.41
  mol_N_P  <- get_nitrogen(P = P) / 14.007
  mol_C_P  <- mol_N_P * 3.7
  carbon_3 <- mol_C_P * 12.011
  unname(carbon_1 + carbon_2 + carbon_3)
}


# ---------------------------------------------------------------------------
# Population functions
# ---------------------------------------------------------------------------

#' Generate a population time-series by back-calculating from harvest
#'
#' Starting from a known harvest number of individuals and a fixed daily
#' mortality rate, integrates the population ODE backwards in time using
#' Euler's method, returning a vector of population size at each time-step
#' in forward chronological order.
#'
#' @param harvest_n Numeric. Number of individuals at harvest (end of grow-out).
#' @param mort      Numeric. Daily mortality rate (proportional mortality, use negative for gains).
#' @param times     Named numeric vector with elements:
#'   \describe{
#'     \item{t_start}{Integer start time-step index.}
#'     \item{t_end}{Integer end time-step index.}
#'     \item{dt}{Time-step size (days).}
#'   }
#'
#' @return Numeric vector of population size at each time-step, from
#'   stocking to harvest.
generate_pop <- function(harvest_n, mort, times) {

  ts <- seq(times['t_start'], times['t_end'], by = times['dt'])   # Integration times

  # Initial condition and vectors initialization
  N_pop <- rep(0, length(ts))                              # Initialize vector N_pop
  N_pop[1] <- harvest_n                                    # Impose harvest condition

  # for cycle that solves population ODE with Euler method
  for (t in 2:length(ts)){
    dN <- unname(mort * N_pop[t - 1])                       # Individuals increment
    N_pop[t] <- N_pop[t - 1] + dN                          # Individuals at time t+1
  }
  return(rev(N_pop))
}


# ---------------------------------------------------------------------------
# Bioenergetic sub-functions
# ---------------------------------------------------------------------------

#' Calculate the temperature-dependent relative feeding rate
#'
#' Implements the H(Tw) function from Baldan et al. (2018), which scales
#' potential ingestion between zero (at thermal limits) and one (at the
#' optimal temperature).
#'
#' @param water_temp   Numeric (scalar or vector). Water temperature (°C).
#' @param species_params Named numeric vector of species parameters. Must
#'   contain \code{betac}, \code{Toa}, and \code{Tma} (see file header for
#'   definitions).
#'
#' @return Numeric. Relative feeding rate in [0, 1].
feeding_rate <- function(water_temp, species_params) {
  exp(species_params['betac'] * (water_temp - species_params['Toa'])) *
    ((species_params['Tma'] - water_temp) / (species_params['Tma'] - species_params['Toa']))^
    (species_params['betac'] * (species_params['Tma'] - species_params['Toa']))
}


#' Calculate the amount of food provided at a time-step
#'
#' If water temperature is within the species' viable feeding range the
#' ration is set to the potential ingestion scaled by a stochastic
#' overfeeding factor; outside that range a baseline minimum ration is
#' returned instead.
#'
#' @param species_params Named numeric vector of species parameters. Must
#'   contain \code{Taa}, \code{Tma}, \code{overFmean}, and \code{overFdelta}
#'   (see file header for definitions).
#' @param water_temp   Numeric. Current water temperature (°C).
#' @param ing_pot      Numeric. Potential ingestion at the current weight and
#'   temperature (g/day).
#' @param ing_pot_min  Numeric. Minimum (baseline) ingestion rate used outside
#'   the viable temperature range (g/day).
#'
#' @return Numeric. Food ration provided (g/day).
food_prov_rate <- function(species_params, water_temp, ing_pot, ing_pot_min) {
  if (water_temp > species_params['Taa'] & water_temp < species_params['Tma']) {
    ing_pot * (1 + rnorm(1, species_params['overFmean'], species_params['overFdelta']))
  } else {
    ing_pot_min
  }
}


#' Apportion provided and ingested feed into waste and assimilated fractions
#'
#' For a single macronutrient component, partitions the provided ration into
#' uneaten feed, faecal excretion, and assimilated material.
#'
#' @param provided           Numeric. Total feed provided (g/day).
#' @param ingested           Numeric. Total feed ingested (g/day).
#' @param ingred_proportion  Numeric. Proportion of this ingredient in the feed
#'   (0–1).
#' @param ingred_macro              Numeric. Energy or mass content per gram of this
#'   macronutrient.
#' @param digestibility Numeric. Digestibility coefficient for this macronutrient
#'   (0–1).
#'
#' @return Named numeric vector with three elements:
#'   \describe{
#'     \item{uneaten}{Mass of uneaten macronutrient (g/day).}
#'     \item{assimilated}{Mass of assimilated macronutrient (g/day).}
#'     \item{excreted}{Mass of faecally excreted macronutrient (g/day).}
#'   }
apportion_feed <- function(provided, ingested, ingred_proportion, ingred_macro, digestibility) {
  provided_g    <- provided * ingred_proportion * ingred_macro
  ingested_g    <- ingested * ingred_proportion * ingred_macro
  assimilated_g <- ingested_g * digestibility

  c(
    uneaten      = sum(provided_g - ingested_g, na.rm = T),
    assimilated  = sum(assimilated_g, na.rm = T),
    excreted     = sum(ingested_g - assimilated_g, na.rm = T)
  )
}


# ---------------------------------------------------------------------------
# Individual fish growth model
# ---------------------------------------------------------------------------

#' Simulate bioenergetic growth of a single individual over time
#'
#' Integrates the full bioenergetic model (based on Baldan et al. 2018) for
#' one individual fish over the specified time period, returning a matrix of
#' state variables at each time-step. Feed composition is supplied via
#' \code{feed_params} and all species biology and population variability
#' parameters are in \code{species_params}.
#'
#' @param species_params Named numeric vector of combined species and population
#'   parameters (see file header for full definitions). Key elements used here
#'   include \code{m}, \code{n}, \code{alpha}, \code{pk}, \code{k0},
#'   \code{omega}, \code{epsO2}, \code{epsprot}, \code{epslip},
#'   \code{epscarb}, \code{a}, \code{k}, \code{eff}, \code{Taa}, \code{Tma},
#'   \code{overFmean}, and \code{overFdelta}.
#' @param water_temp     Numeric vector. Water temperature (°C) at each
#'   time-step. Length must equal the number of days in the simulation.
#' @param feed_params    Named list with one element per macronutrient
#'   (\code{Proteins}, \code{Lipids}, \code{Carbohydrates}). Each element is
#'   itself a list with fields \code{ingred_proportion}, \code{ingred_macro}, and
#'   \code{digest}.
#' @param times          Named numeric vector with elements \code{t_start},
#'   \code{t_end}, and \code{dt} (time-step size in days).
#' @param init_weight    Numeric. Initial individual weight (g).
#' @param ingmax         Numeric. Maximum mass-specific ingestion coefficient
#'   (g g\eqn{^{-m}} day\eqn{^{-1}}).
#' @param output_vars    Character vector. Names of variables to include in the
#'   returned matrix. Any subset of the 22 base variables (\code{weight},
#'   \code{dw}, \code{water_temp}, \code{T_response}, \code{P_excr},
#'   \code{L_excr}, \code{C_excr}, \code{P_uneat}, \code{L_uneat},
#'   \code{C_uneat}, \code{food_prov}, \code{food_enc}, \code{rel_feeding},
#'   \code{ing_pot}, \code{ing_act}, \code{E_assim}, \code{E_somat},
#'   \code{anab}, \code{catab}, \code{O2}, \code{NH4}) and/or the 9 derived
#'   aggregate variables (\code{total_excr}, \code{total_uneat}, \code{metab},
#'   \code{nitrogen_excr}, \code{nitrogen_uneat}, \code{carbon_excr},
#'   \code{carbon_uneat}, \code{total_carbon}, \code{total_nitrogen}). Defaults
#'   to all 31 variables.
#'
#' @return Numeric matrix with one row per day. The \code{days} column is always
#'   present; remaining columns are those specified in \code{output_vars}.
fish_growth <- function(
  species_params,
  water_temp,
  feed_params,
  times,
  init_weight,
  ingmax,
  output_vars = c("weight", "dw", "water_temp", "T_response", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "food_prov", "food_enc", "rel_feeding", "ing_pot", "ing_act", "E_assim", "E_somat", "anab", "catab", "O2", "NH4", "total_excr", "total_uneat", "metab", "nitrogen_excr", "nitrogen_uneat", "carbon_excr", "carbon_uneat", "total_carbon", "total_nitrogen")
) {
  # Pre-calculate array sizes
  n_days <- length(times['t_start']:times['t_end'])

  # Always allocate all 22 base columns — the computation loop requires all of
  # them regardless of output_vars. Aggregate columns are appended after the loop.
  base_cols <- c('days', 'weight', 'dw', 'water_temp', 'T_response', 'P_excr',
                 'L_excr', 'C_excr', 'P_uneat', 'L_uneat', 'C_uneat', 'food_prov',
                 'food_enc', 'rel_feeding', 'ing_pot', 'ing_act', 'E_assim',
                 'E_somat', 'anab', 'catab', 'O2', 'NH4')
  result <- matrix(NA, nrow = n_days, ncol = length(base_cols))
  colnames(result) <- base_cols

  # Initialize first values
  result[, 'days']       <- (times['t_start']:times['t_end']) * times['dt']
  result[1, 'weight']    <- init_weight
  result[, 'water_temp'] <- water_temp

  # Main calculation loop
  for (i in 1:(n_days - 1)) {
    # Temperature response and ingestion calculations
    result[i, 'rel_feeding'] <- feeding_rate(result[i, 'water_temp'], species_params)
    result[i, 'ing_pot']     <- ingmax * (result[i, 'weight']^species_params['m']) * result[i, 'rel_feeding']

    # Food provision and ingestion
    result[i, 'food_prov'] <- food_prov_rate(
      species_params = species_params,
      water_temp     = result[i, 'water_temp'],
      ing_pot        = result[i, 'ing_pot'],
      ing_pot_min    = ingmax * (result[i, 'weight']^species_params['m']) *
                         feeding_rate(species_params['Taa'], species_params)
    )
    result[i, 'food_enc'] <- species_params['eff'] * result[i, 'food_prov']
    result[i, 'ing_act']  <- min(result[i, 'food_enc'], result[i, 'ing_pot'])

    # Energy content of somatic tissue at current weight
    result[i, 'E_somat'] <- species_params['a'] * result[i, 'weight']^species_params['k']

    # Process feed components
    app_carbs    <- apportion_feed(
      result[i, 'food_prov'],
      result[i, 'ing_act'],
      feed_params[['Carbohydrates']]$proportion,
      feed_params[['Carbohydrates']]$macro,
      feed_params[['Carbohydrates']]$digest
    )
    app_lipids   <- apportion_feed(
      result[i, 'food_prov'], result[i, 'ing_act'],
      feed_params[['Lipids']]$proportion,
      feed_params[['Lipids']]$macro,
      feed_params[['Lipids']]$digest
    )
    app_proteins <- apportion_feed(
      result[i, 'food_prov'], result[i, 'ing_act'],
      feed_params[['Proteins']]$proportion,
      feed_params[['Proteins']]$macro,
      feed_params[['Proteins']]$digest
    )

    # Excretion and uneaten-feed waste
    result[i, c('C_excr', 'L_excr', 'P_excr')]   <- c(app_carbs['excreted'],  app_lipids['excreted'],  app_proteins['excreted'])
    result[i, c('C_uneat', 'L_uneat', 'P_uneat')] <- c(app_carbs['uneaten'],   app_lipids['uneaten'],   app_proteins['uneaten'])

    # Energy assimilated from feed
    result[i, 'E_assim'] <-
      app_carbs['assimilated']    * species_params['epscarb'] +
      app_lipids['assimilated']   * species_params['epslip']  +
      app_proteins['assimilated'] * species_params['epsprot']

    # Temperature-dependent metabolism
    result[i, 'T_response'] <- exp(species_params['pk'] * result[i, 'water_temp'])
    result[i, 'anab']       <- result[i, 'E_assim'] * (1 - species_params['alpha'])
    result[i, 'catab']      <- species_params['epsO2'] * species_params['k0'] *
      result[i, 'T_response'] * (result[i, 'weight']^species_params['n']) *
      species_params['omega']

    # Oxygen consumption and ammonium excretion
    result[i, 'O2']  <- result[i, 'catab'] / species_params['epsO2']
    result[i, 'NH4'] <- result[i, 'O2'] * 0.06

    # Weight increment (Euler integration)
    result[i, 'dw']           <- (result[i, 'anab'] - result[i, 'catab']) / result[i, 'E_somat']
    result[i + 1, 'weight']   <- result[i, 'weight'] + result[i, 'dw'] * times['dt']
  }

  # Always compute all aggregate variables — they all source from base columns
  # which are guaranteed to be present. This ensures dependency chains (e.g.
  # total_carbon requires carbon_excr) are always satisfied regardless of what
  # the user puts in output_vars. Subset to output_vars only at return.
  nitrogen_excr  <- get_nitrogen(result[, 'P_excr'])
  nitrogen_uneat <- get_nitrogen(result[, 'P_uneat'])
  carbon_excr    <- get_carbon(P = result[, 'P_excr'],  L = result[, 'L_excr'],  C = result[, 'C_excr'])
  carbon_uneat   <- get_carbon(P = result[, 'P_uneat'], L = result[, 'L_uneat'], C = result[, 'C_uneat'])

  agg_mat <- cbind(
    total_excr     = result[, 'P_excr']  + result[, 'L_excr']  + result[, 'C_excr'],
    total_uneat    = result[, 'P_uneat'] + result[, 'L_uneat'] + result[, 'C_uneat'],
    metab          = result[, 'anab'] - result[, 'catab'],
    nitrogen_excr  = nitrogen_excr,
    nitrogen_uneat = nitrogen_uneat,
    carbon_excr    = carbon_excr,
    carbon_uneat   = carbon_uneat,
    total_carbon   = carbon_excr   + carbon_uneat,
    total_nitrogen = nitrogen_excr + nitrogen_uneat
  )

  result <- cbind(result, agg_mat)

  # Subset to days + the user-requested output_vars
  keep_cols <- c('days', base::intersect(output_vars, colnames(result)))
  return(result[, keep_cols, drop = FALSE])
}


# ---------------------------------------------------------------------------
# Farm-level growth models
# ---------------------------------------------------------------------------

#' Simulate all individuals in a farm and return full per-individual trajectories
#'
#' The core farm simulation function. Runs \code{fish_growth()} for \code{MC_pop}
#' individuals in parallel (via \pkg{furrr}) and collects each individual's full
#' time-series into a list of \code{MC_pop × n_days} matrices. Derived aggregate
#' and elemental variables are computed at the farm level (more efficient than
#' per-individual calculation). No consolidation to farm-level summaries or
#' population scaling is applied — use \code{consolidate_farm()} and
#' \code{scale_farm_population()} for those steps, or the \code{farm_growth()}
#' wrapper which calls all three in sequence.
#'
#' When \code{use_MC_population = TRUE} (default), individual initial weights and
#' maximum ingestion rates are drawn from normal distributions parameterised by
#' \code{meanW}/\code{deltaW} and \code{meanImax}/\code{deltaImax} in
#' \code{species_params}. When \code{use_MC_population = FALSE}, all individuals
#' are identical (variability parameters are zeroed) and \code{MC_pop} is forced
#' to 1.
#'
#' @param species_params    Named numeric vector of combined species and
#'   population parameters (see file header). Must include \code{meanW},
#'   \code{deltaW}, \code{meanImax}, and \code{deltaImax} in addition to all
#'   bioenergetic parameters required by \code{fish_growth()}.
#' @param feed_params       Named list of macronutrient parameters (see
#'   \code{fish_growth()}).
#' @param water_temp        Numeric vector. Water temperature (°C) at each
#'   time-step.
#' @param times             Named numeric vector with \code{t_start},
#'   \code{t_end}, and \code{dt}.
#' @param use_MC_population Logical. If \code{TRUE} (default), sample individual
#'   weights and ingestion rates from normal distributions to represent
#'   population variability. If \code{FALSE}, all individuals are identical.
#' @param MC_pop            Integer. Number of individuals to simulate. Ignored
#'   (forced to 1) when \code{use_MC_population = FALSE}.
#' @param output_vars       Character vector. Names of variables to return (see
#'   \code{fish_growth()} for the full list of available base and aggregate
#'   variables). Defaults to all 31 variables. The \code{days} variable is always
#'   retained internally for use by \code{consolidate_farm()}.
#'
#' @return Named list of matrices, one per variable in \code{output_vars} plus
#'   \code{days}, each of dimensions \code{MC_pop × n_days} (individuals as
#'   rows, time-steps as columns).
farm_growth_full <- function(
  species_params, 
  feed_params, 
  water_temp, 
  times,
  use_MC_population = TRUE, 
  MC_pop, 
  output_vars = c("weight", "dw", "water_temp", "T_response", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "food_prov", "food_enc", "rel_feeding", "ing_pot", "ing_act", "E_assim", "E_somat", "anab", "catab", "O2", "NH4", "total_excr", "total_uneat", "metab", "nitrogen_excr", "nitrogen_uneat", "carbon_excr", "carbon_uneat", "total_carbon", "total_nitrogen")
) {

  # When uniform population is requested, zero out variability and use one individual
  if (!use_MC_population) {
    species_params['deltaW']    <- 0
    species_params['deltaImax'] <- 0
    MC_pop <- 1
  }

  # Generate all random values upfront
  init_weights <- rnorm(MC_pop, mean = species_params['meanW'],    sd = species_params['deltaW'])
  ingmaxes     <- rnorm(MC_pop, mean = species_params['meanImax'], sd = species_params['deltaImax'])

  # Run parallel simulation for individuals
  mc_results <- future_map2(init_weights, ingmaxes, function(init_w, ing_m) {
    fish_growth(
      species_params = species_params,
      water_temp     = water_temp,
      feed_params    = feed_params,
      times          = times,
      init_weight    = init_w,
      ingmax         = ing_m,
      output_vars    = output_vars
    )
  },
  .options = furrr_options(seed = TRUE)
  )

  stat_names <- colnames(mc_results[[1]])
  # Stack all individuals — individuals are rows, time-steps are columns
  all_results <- lapply(1:length(stat_names), function(col_idx) {
    t(sapply(mc_results, function(mat_idx) mat_idx[, col_idx]))
  }) %>%
    setNames(stat_names)

  # Subset to only the requested output_vars (days always kept for consolidate_farm())
  return(all_results)
}


#' Consolidate individual fish results to farm-level means and standard deviations
#'
#' Takes the list-of-matrices output from \code{farm_growth_full()} and reduces
#' each variable across individuals, producing a per-time-step mean and standard
#' deviation. The \code{days} variable is extracted and returned separately as it
#' is the same for all individuals and does not require summarising.
#'
#' @param full_results Named list of \code{MC_pop × n_days} matrices as returned
#'   by \code{farm_growth_full()}.
#'
#' @return Named list with two elements:
#'   \describe{
#'     \item{days}{Numeric vector of time values (days).}
#'     \item{stats}{Named list of \code{n_days × 2} matrices, one per variable
#'       (excluding \code{days}). Column 1 is the cross-individual mean;
#'       column 2 is the cross-individual standard deviation.}
#'   }
consolidate_farm <- function(full_results) {

  days      <- full_results[["days"]][1, ]   # identical for all individuals; take first row
  var_names <- names(full_results)[names(full_results) != "days"]

  stats <- lapply(var_names, function(nm) {
    cbind(
      colMeans(full_results[[nm]]),
      matrixStats::colSds(full_results[[nm]])
    ) %>%
      as.matrix() %>% unname()
  }) %>%
    setNames(var_names)

  list(days = days, stats = stats)
}


#' Scale population-level farm outputs by the stocking/mortality trajectory
#'
#' Takes the output of \code{consolidate_farm()} and multiplies a defined set of
#' population-total variables (waste fluxes, elemental outputs, O2, NH4, etc.)
#' by \code{N_pop} at each time-step. The standard deviation column for each
#' scaled variable is converted from absolute to relative (CV) before scaling,
#' then back to absolute afterwards. The result mirrors the structure returned by
#' \code{consolidate_farm()}.
#'
#' @param consolidated Named list as returned by \code{consolidate_farm()}, with
#'   elements \code{days} (numeric vector) and \code{stats} (named list of
#'   \code{n_days × 2} matrices).
#' @param N_pop        Numeric vector. Number of individuals in the farm at each
#'   time-step (from \code{generate_pop()}). Must be at least as long as the
#'   number of time-steps.
#' @param scaled_vars  Character vector. Names of variables in \code{stats} to
#'   scale by \code{N_pop} and return. Only variables present in both
#'   \code{scaled_vars} and \code{consolidated$stats} are processed. Defaults to
#'   a selection of population-scalable variables: \code{weight}, \code{P_excr}, 
#'   \code{L_excr}, \code{C_excr}, \code{P_uneat}, \code{L_uneat}, \code{C_uneat},
#'   \code{total_excr}, \code{total_uneat}, \code{NH4}, \code{food_prov},
#'   \code{nitrogen_excr}, \code{nitrogen_uneat}, \code{carbon_excr},
#'   \code{carbon_uneat}, \code{total_carbon}, \code{total_nitrogen}.
#'
#' @return Named list with two elements:
#'   \describe{
#'     \item{days}{Numeric vector of time values (days).}
#'     \item{stats}{Named list of \code{n_days × 2} matrices, one per variable
#'       in \code{scaled_vars} (intersected with available variables), with names
#'       suffixed by \code{_scaled}. Column 1 is the population-scaled mean;
#'       column 2 is the population-scaled standard deviation.}
#'   }
scale_farm_population <- function(
  consolidated, 
  N_pop, 
  scaled_vars = c("weight", "P_excr", "L_excr", "C_excr", "P_uneat", 
  "L_uneat", "C_uneat", "total_excr", "total_uneat", "NH4", 
  "food_prov", "nitrogen_excr", "nitrogen_uneat", "carbon_excr", "carbon_uneat", 
  "total_carbon", "total_nitrogen")
) {

  days  <- consolidated[["days"]]
  stats <- consolidated[["stats"]]

  # Variables whose values represent per-individual means and must be
  # multiplied by N_pop to obtain farm-level totals
  for (stat_nm in scaled_vars) {
    # Convert SD from absolute to relative (CV), scale mean, then back to absolute
    stats[[stat_nm]][, 2] <- stats[[stat_nm]][, 2] / stats[[stat_nm]][, 1]
    stats[[stat_nm]][, 1] <- stats[[stat_nm]][, 1] * N_pop[1:length(days)]
    stats[[stat_nm]][, 2] <- stats[[stat_nm]][, 1] * stats[[stat_nm]][, 2]
  }

  # Subset to only the scaled variables, rename with _scaled suffix,
  # and return in consolidate_farm() structure
  stats <- stats[base::intersect(scaled_vars, names(stats))]
  names(stats) <- paste0(names(stats), "_scaled")

  list(days = days, stats = stats)
}


#' Simulate farm-level growth with optional Monte Carlo population variability
#'
#' Convenience wrapper that calls \code{farm_growth_full()},
#' \code{consolidate_farm()}, and \code{scale_farm_population()} in sequence,
#' returning consolidated, population-scaled farm output. When
#' \code{use_MC_population = TRUE} (default), individual initial weights and
#' maximum ingestion rates are drawn from normal distributions; when
#' \code{FALSE}, all individuals are identical and \code{MC_pop} is forced to 1.
#'
#' @param species_params    Named numeric vector of combined species and
#'   population parameters (see file header). Must include \code{meanW},
#'   \code{deltaW}, \code{meanImax}, and \code{deltaImax} in addition to all
#'   bioenergetic parameters required by \code{fish_growth()}.
#' @param feed_params       Named list of macronutrient parameters (see
#'   \code{fish_growth()}).
#' @param water_temp        Numeric vector. Water temperature (°C) at each
#'   time-step.
#' @param times             Named numeric vector with \code{t_start},
#'   \code{t_end}, and \code{dt}.
#' @param N_pop             Numeric vector. Number of individuals in the farm
#'   at each time-step (from \code{generate_pop()}).
#' @param use_MC_population Logical. If \code{TRUE} (default), sample
#'   individual weights and ingestion rates from normal distributions to
#'   represent population variability. If \code{FALSE}, all individuals are
#'   identical.
#' @param MC_pop            Integer. Number of Monte Carlo individuals to
#'   simulate. Ignored (forced to 1) when \code{use_MC_population = FALSE}.
#' @param output_vars       Character vector or \code{NULL}. Names of variables
#'   to simulate and return (passed to \code{farm_growth_full()}; see
#'   \code{fish_growth()} for the full variable list). \code{NULL} (default)
#'   returns all 31 variables.
#' @param scaled_vars       Character vector. Names of variables to scale by
#'   \code{N_pop} (passed to \code{scale_farm_population()}). Defaults to all
#'   population-scalable variables (see \code{scale_farm_population()}).
#'
#' @return Named list with two elements:
#'   \describe{
#'     \item{days}{Numeric vector of time values (days).}
#'     \item{stats}{Named list of \code{n_days × 2} matrices combining the
#'       outputs of \code{consolidate_farm()} (per-individual means/SDs for
#'       variables in \code{output_vars}) and \code{scale_farm_population()}
#'       (population-scaled means/SDs for variables in \code{scaled_vars},
#'       with names suffixed by \code{_scaled}).}
#'   }
farm_growth <- function(
  species_params,
  feed_params,
  water_temp,
  times,
  N_pop,
  use_MC_population = TRUE,
  MC_pop,
  output_vars = NULL,
  scaled_vars = c("weight", "dw", "P_excr", "L_excr", "C_excr", "P_uneat", "L_uneat", "C_uneat", "ing_act", "total_excr", "total_uneat", "O2", "NH4", "food_prov", "nitrogen_excr", "nitrogen_uneat", "carbon_excr", "carbon_uneat", "total_carbon", "total_nitrogen")
) {
  # NULL output_vars defers to farm_growth_full()'s own default (all variables)
  if (is.null(output_vars)) {
    full_results <- farm_growth_full(species_params, feed_params, water_temp, times, use_MC_population, MC_pop)
  } else {
    full_results <- farm_growth_full(species_params, feed_params, water_temp, times, use_MC_population, MC_pop, output_vars = output_vars)
  }
  consolidated  <- consolidate_farm(full_results)
  scaled        <- scale_farm_population(consolidated, N_pop, scaled_vars = scaled_vars)

  list(
    days  = consolidated[["days"]],
    stats = c(consolidated[["stats"]], scaled[["stats"]])
  )
}


# ---------------------------------------------------------------------------
# Utility / post-processing
# ---------------------------------------------------------------------------

#' Subset and reshape full farm results to a tidy long-format data frame
#'
#' Takes the list-of-matrices output from \code{farm_growth_full()} and
#' converts a selected subset of variables into a single tidy (long-format)
#' data frame suitable for plotting and analysis with \pkg{ggplot2} or
#' \pkg{dplyr}.
#'
#' @param full_results Named list of matrices as returned by
#'   \code{farm_growth_full()}. Each matrix is \code{MC_pop} × \code{n_days}.
#' @param out_vars Character vector. Names of variables to retain. Must be
#'   valid names present in \code{full_results}. Available variables include
#'   the 22 base outputs from \code{fish_growth()} plus the 9 derived aggregate
#'   variables: \code{total_excr}, \code{total_uneat}, \code{metab},
#'   \code{nitrogen_excr}, \code{nitrogen_uneat}, \code{carbon_excr},
#'   \code{carbon_uneat}, \code{total_carbon}, \code{total_nitrogen}.
#'
#' @return A tidy long-format data frame with columns:
#'   \describe{
#'     \item{fish}{Integer index identifying each simulated individual.}
#'     \item{prod_t}{Integer time-step index.}
#'     \item{value}{Numeric value of the output variable.}
#'     \item{measure}{Factor indicating the output variable name.}
#'   }
decomposed_to_tidy <- function(
  full_results, 
  output_vars = NULL
) {

  selected_vars <- if (is.null(output_vars)) {
    names(full_results)[names(full_results) != "days"]
  } else {
    output_vars
  }

  subset_results <- full_results[selected_vars]

  map_dfr(1:length(subset_results), function(col_idx) {
    res <- reshape2::melt(subset_results[[col_idx]]) %>%
      mutate(measure = as.factor(selected_vars[col_idx]))
    colnames(res) <- c("fish", "prod_t", "value", "measure")
    res
  })
}
