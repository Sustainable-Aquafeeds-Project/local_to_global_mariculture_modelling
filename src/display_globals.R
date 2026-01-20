suppressPackageStartupMessages(suppressWarnings({
  library(tidyverse)
  library(RColorBrewer)
  library(colorspace)
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

orig_feeds <- c("marine_dominant_biomar", "animal_inclusive_biomar", "novel_inclusive_biomar")
long_feeds <- c("marine_dominant_biomar", "animal_dominant_biomar", "plant_dominant_biomar")
short_feeds <- c("MD", "AD", "PD")

# base_cols <- brewer.pal(3, "Dark2")
# "#1B9E77" "#D95F02" "#7570B3"

feed_pal <- c(
  "PD" = "#1B9E77",
  "MD" = "#7570B3",
  "AD" = "#D95F02"#,
  # "plant_dominant_biomar" = "#1B9E77",
  # "marine_dominant_biomar" = "#7570B3",
  # "animal_dominant_biomar" = "#D95F02"
)
cohort_pal <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", "all" = "black")

macro_pal <- c(
  "protein" = "#c73535ff", 
  "P" = "#c73535ff", 
  "carbohydrates" = "#4682B4", 
  "C" = "#4682B4", 
  "lipids" = "#FFB90F", 
  "L" = "#FFB90F"
)

# lighten("steelblue", amount = c(0.35, 0, -0.35))
# "#81B4E7" "#4682B4" "#054E77"

level_pal <- c(
  "min" = "#81B4E7", "Minimum" = "#81B4E7",
  "mean" = "#4682B4", "Mean" = "#4682B4",
  "max" = "#054E77", "Maximum" = "#054E77"
)
level_lines <- c(
  "min" = "dotted", "Minimum" = "dotted",
  "mean" = "solid", "Mean" = "solid",
  "max" = "dashed", "Maximum" = "dashed"
)
level_shapes <- c("Mean" = 21, "Minimum" = 2, "Maximum" = 4)

