library(devtools)
library(qs)
library(arrow)
library(stringr)

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

no_margins <- function() {
  theme(legend.position = "none", 
        plot.margin = margin(0, 0, 0, 0),
        axis.title = element_blank())
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

prettyplot <- function() {
  theme_classic() +
    theme(legend.position = "none",
          text = element_text(family = "serif", size = 12, colour = "black"),
          axis.title.x = element_text(vjust = 0.5),
          axis.title.y = element_text(hjust = 0.5))
}
