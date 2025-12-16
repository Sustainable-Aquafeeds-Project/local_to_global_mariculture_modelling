suppressPackageStartupMessages(suppressWarnings({
  library(tidyverse)
}))

# This function is designed to work with a DPI of 300
prettyplot_300 <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "serif", size = 12, colour = "black"),
      legend.position = "none",
      axis.title.y = element_text(vjust = 0.5),
      axis.title.x = element_text(vjust = 0.5),
      legend.title = element_blank()
    )
}

# This function is designed to work with a DPI of 600
prettyplot_600 <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "serif", size = 12, colour = "black"),
      legend.position = "none",
      axis.title.y = element_text(vjust = 0.5),
      axis.title.x = element_text(vjust = 0.5),
      legend.title = element_blank()
    )
}

long_feeds <- c("marine_dominant_biomar", "plant_dominant", "novel_inclusive")
short_feeds <- c("MD", "PD", "NI")

feed_pal <- c(
  "past" = "#E41A1C", 
  "reference" = "#377EB8", 
  "future" = "#4DAF4A",
  "plant_dominant" = "#4DAF4A",
  "marine_dominant" = "#377EB8",
  "novel_inclusive" = "#E41A1C",
  "PD" = "#4DAF4A",
  "MD" = "#377EB8",
  "NI" = "#E41A1C"
)
cohort_pal <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", "all" = "black")
macro_pal <- c(
  "protein" = "darkred", 
  "P" = "darkred", 
  "carbohydrates" = "steelblue", 
  "C" = "steelblue", 
  "lipids" = "darkgoldenrod1", 
  "L" = "darkgoldenrod1"
)
