# nolint start

# Setup -----------------------------------------------------------------------------------------------------------
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
conflicts_prefer(dplyr::filter(), dplyr::select(), .quiet = T)

here("00_model_functions.R") %>% source()
here("00_dirs.R") %>% source()

ingred_inputs_file <- file.path(input_feed_profile_path, "all_ingredients.csv")
feed_input_file <- file.path(input_feed_profile_path, "all_feeds.csv")
feed_params_file <- file.path(output_species_data_path, "feed_params.qs")

# Ingredients -----------------------------------------------------------------------------------------------------
# Warning: there is no check here to ensure that protein+lipid+ash+carb = 1 in the incoming data. 
# You must make sure that the proportions add up to 100% manually.
ingreds <- ingred_inputs_file %>% 
  read.csv() %>% 
  mutate(ingredient = as.factor(ingredient))
ingred_nms <- levels(ingreds$ingredient)

# Feeds -----------------------------------------------------------------------------------------------------------
feed_inputs <- feed_input_file %>% 
  read.csv() %>% 
  mutate(feed = as.factor(feed),
         ingredient = as.factor(ingredient)) %>% 
  merge(ingreds, by = "ingredient", all = T)
feed_types <- levels(feed_inputs$feed)

feed_params <- purrr::map(feed_types, function(ft) {
  df <- feed_inputs %>% 
    filter(feed == ft) 
  list(
    Proteins = df %>% 
      select(ingredient, proportion, contains("protein"), -contains("feed")) %>%
      rename(macro = protein, digest = protein_digestibility),
    Carbohydrates = df %>% 
      select(ingredient, proportion, contains("carb"), -contains("feed")) %>%
      rename(macro = carb, digest = carb_digestibility),
    Lipids = df %>% 
      select(ingredient, proportion, contains("lipid"), -contains("feed")) %>%
      rename(macro = lipid, digest = lipid_digestibility)
  )
}) %>% 
  setNames(feed_types)

qsave(feed_params, feed_params_file)
