suppressPackageStartupMessages(suppressWarnings({
library(devtools)
library(qs)
library(dplyr)
library(magrittr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(purrr)
}))

# Define areas -------------------------------------------------
# For drawing area boxes on the big map
inset_boxes <- list(   # lonmin, lonmax, latmin, latmax
  CAN1 = c(-132, -122, 47.75, 54.1),
  CAN2 = c(-69, -54.75, 43.25, 48.25),
  EUR = c(-26, 30, 51, 72),
  CHI = c(-77.5, -62.5, -56, -25),
  AUS = c(144, 149.5, -44, -39.75)
)

# For the corresponding patchwork map lims, basically the same, just a little refined
inset_boxes_sm <- list(   # lonmin, lonmax, latmin, latmax
  CAN1 = c(-130.5, -123, 48.25, 54),
  CAN2 = c(-69, -55, 43.5, 48),
  EUR = c(-23, 28.5, 52, 71),
  CHI = c(-77, -62, -55.25, -27.5),
  AUS = c(144.5, 148.5, -43.75, -40.75)
) %>% 
  map(function(bx) {
    list(
      xlims = bx[1:2],
      ylims = bx[3:4],
      labx = bx[1],
      laby = bx[4]
    )
})

# Where should the inset labels be positioned in the big map
labels_spec_robinson <- c(
  CAN1 = "bottom_left_outside", 
  CAN2 = "bottom_right_outside", 
  EUR = "bottom_left_outside", 
  CHI = "top_left_outside", 
  AUS = "bottom_left_outside"
)

labels_offset_robinson <- c(
  CAN1 = 1, 
  CAN2 = 2, 
  EUR = 1.75, 
  CHI = 1.75, 
  AUS = 3.5
)

labels_spec_mercator <- list(
  CAN1 = c(l = "A", h = 0, v = 0), 
  CAN2 = c(l = "B", h = 0.75, v = 0), 
  EUR = c(l = "C", h = -0.75, v = 0.5), 
  CHI = c(l = "D", h = 0, v = 0), 
  AUS = c(l = "E", h = 0, v = 0)
)

# The Robinson projection --------------------------------------------------------------------------------------------------------------
worldmap_robinson <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_transform(crs = "+proj=robin")
graticules_robinson <- st_graticule(worldmap_robinson, lon = seq(-180, 180, 30), lat = seq(-90, 90, 30)) 

box_data <- create_boxes(
  box_list = inset_boxes, # see function create_boxes above
  label_positions = unname(labels_spec_robinson),
  offset_deg = unname(labels_offset_robinson)
  )

boxes_robinson <- st_transform(box_data$boxes, crs = "+proj=robin")
labels_robinson <- st_transform(box_data$labels, crs = "+proj=robin")

p_bigmap_robinson <- ggplot() +
  geom_sf(data = graticules_robinson, color = "gray80", size = 0.3) +
  geom_sf(data = worldmap_robinson, fill = "white", color = "dimgray") +
  coord_sf() +
  theme_void()

p_bigmap_robinson_boxes <- ggplot() +
  geom_sf(data = graticules_robinson, color = "gray80", size = 0.3) +
  geom_sf(data = worldmap_robinson, fill = "white", color = "dimgray") +
  coord_sf() +
  geom_sf(data = boxes_robinson, fill = NA, color = "darkred", size = 1.65) +
  geom_sf_text(data = labels_robinson, aes(label = letter), 
               color = "darkred", size = 5, fontface = "bold", 
               hjust = 0.5, vjust = 0.5) +
  theme_void()

# The Mercator projection --------------------------------------------------------------------------------------------------------------
worldmap_mercator <- ne_countries(scale = "large", returnclass = "sf")

get_insets_mercator <- function(map) {
  map(inset_boxes_sm, function(specs) {
    map + 
      coord_sf(
        xlim = specs[["xlims"]], 
        ylim = specs[["ylims"]]
      )
  })
}

p_bigmap_mercator <- ggplot() +
  geom_sf(data = worldmap_mercator, fill = "white", color = "dimgray") +
  coord_sf() +
  labs(y = "Latitude", x = "Longitude")

p_mercator_insets <- get_insets_mercator(p_bigmap_mercator)

patchwork_mercator <- function(bigmap) {
  p_insets <- get_insets_mercator(bigmap)
  CAN1 <- p_insets[["CAN1"]] + scale_x_continuous(breaks = seq(-130, -122, 2))
  CAN2 <- p_insets[["CAN2"]] + scale_x_continuous(breaks = seq(-72.5, -50, 5))
  EUR <- p_insets[["EUR"]] + scale_x_continuous(breaks = seq(-30, 50, 10))
  AUS <- p_insets[["AUS"]] + scale_x_continuous(breaks = seq(140, 150, 1))
  CHI <- p_insets[["CHI"]] + scale_x_continuous(breaks = seq(-80, -60, 5))

  grid_canada <- plot_grid(
    CAN1, CAN2, 
    ncol = 2, 
    rel_widths =  c(1, 1.3)
  )

  grid_notcanada <- plot_grid(
    plot_grid(
      EUR, AUS,  
      nrow = 2, 
      rel_heights = c(1, 1)
    ), 
    CHI,  
    ncol = 2, 
    rel_widths =  c(1.1, 1)
  )
  
  plot_grid(
    grid_canada, grid_notcanada,
    nrow = 2,
    rel_heights = c(1, 1.5)
  )
}

patchwork_mercator_1 <- function(bigmap) {
  p_insets <- get_insets_mercator(bigmap)
  CAN1 <- p_insets[["CAN1"]] + 
    scale_x_continuous(breaks = seq(-130, -122, 2)) +
    theme(legend.position = "none")
  CAN2 <- p_insets[["CAN2"]] + 
    scale_x_continuous(breaks = seq(-70, -50, 2.5)) +
    theme(legend.position = "top")

  plot_grid(
    CAN1, CAN2, 
    ncol = 2, 
    rel_widths =  c(1, 1.65)
  )
}

patchwork_mercator_2 <- function(bigmap) {
  p_insets <- get_insets_mercator(bigmap)
  EUR <- p_insets[["EUR"]] + 
    scale_x_continuous(breaks = seq(-30, 50, 10)) +
    theme(legend.position = "none")
  AUS <- p_insets[["AUS"]] + 
    scale_x_continuous(breaks = seq(140, 150, 1)) +
    theme(legend.position = "none")
  CHI <- p_insets[["CHI"]] + 
    scale_x_continuous(breaks = seq(-82.5, -60, 5)) +
    theme(legend.position = "none")

  plot_grid(
    plot_grid(
      EUR, AUS,  
      nrow = 2, 
      rel_heights = c(1, 0.9)
    ), 
    CHI,  
    ncol = 2, 
    rel_widths =  c(1.1, 1)
  )
}

# patchwork_mercator(p_bigmap_mercator)
# patchwork_mercator_1(p_bigmap_mercator)
# patchwork_mercator_2(p_bigmap_mercator)

# Other map display functions
remove_map_axes <- function() {
  theme(
    axis.title.x = element_blank(), axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )
}
