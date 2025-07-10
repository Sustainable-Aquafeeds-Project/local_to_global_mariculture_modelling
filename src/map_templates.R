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

# Define areas -------------------------------------------------
# For drawing area boxes on the big map
inset_boxes <- list(   # lonmin, lonmax, latmin, latmax
  CAN1 = c(-132, -122, 47.75, 54.25),
  CAN2 = c(-70, -54, 43, 48.5),
  EUR = c(-26, 30, 51, 72),
  CHI = c(-77.5, -62.5, -56, -25),
  AUS = c(144, 149.5, -44, -39.75)
)

# For the corresponding patchwork map lims
inset_boxes_sm <- list(   # lonmin, lonmax, latmin, latmax
  CAN1 = c(-132, -122, 48.5, 54),
  CAN2 = c(-70, -54, 43, 48.5),
  EUR = c(-23, 28.5, 52, 71),
  CHI = c(-78, -62, -55.5, -27),
  AUS = c(144, 149, -43.75, -40.75)
) %>% 
  map(function(bx) {
    list(
      xlims = bx[1:2],
      ylims = bx[3:4],
      labx = bx[1],
      laby = bx[4]
      )
})

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
  CAN1 = c(h = -1, v = 0.5), 
  CAN2 = c(h = 0.5, v = 0.5), 
  EUR = c(h = -1, v = 0.5), 
  CHI = c(h = 0, v = 0), 
  AUS = c(h = 0, v = 0.5)
)

# The Robinson projection --------------------------------------
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
  geom_sf(data = boxes_robinson, fill = NA, color = "darkred", size = 1.65) +
  geom_sf_text(data = labels_robinson, aes(label = letter), 
               color = "darkred", size = 5, fontface = "bold", 
               hjust = 0.5, vjust = 0.5) +
  theme_void()

# The Mercator projection --------------------------------------
worldmap_mercator <- ne_countries(scale = "large", returnclass = "sf")

get_insets_mercator <- function(map) {
  map2(inset_boxes_sm, labels_spec_mercator, function(specs, labs) {
    map + 
      coord_sf(xlim = specs[["xlims"]], 
              ylim = specs[["ylims"]])  +
      draw_label("A", size = 14, fontface = "bold", 
                x = specs[["labx"]], 
                y = specs[["laby"]], 
                hjust = labs["h"], vjust = labs["h"])
  })
}

p_bigmap_mercator <- ggplot() +
  geom_sf(data = worldmap_mercator, fill = "white", color = "dimgray") +
  coord_sf() +
  theme_void() +
  labs(y = "Latitude", x = "Longitude")

p_insets_mercator <- get_insets_mercator(p_bigmap_mercator)


