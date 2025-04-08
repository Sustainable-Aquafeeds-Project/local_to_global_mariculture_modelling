#Template rasters


base_rast <- rast(res = 0.5)
values(base_rast) <- 1:ncell(base_rast)


gall_rast <- rast(ext(-20037508, 20037508, -10000000, 10000000), 
                  resolution = 10000, 
                  crs = "ESRI:54016")

values(gall_rast) <- 1:ncell(gall_rast)

#bounding boxes
east_canada_bbox <- c(xmin = -75, xmax = -45, ymin = 39, ymax = 54)
west_canada_bbox <- c(xmin =-140, xmax = -115, ymin = 45, ymax = 57)
chile_bbox <- c(xmin = -82, xmax = -63, ymin = -57, ymax = -18)
europe_iceland_bbox <- c(xmin = -29, xmax = 48, ymin = 46, ymax = 73)
australia_bbox <- c(xmin = 143, xmax = 149, ymin = -44, ymax = -40)