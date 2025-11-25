# R/calculate_lc_density.R

#' Calculates the focal density of specific Land Cover classes.
#'
#' This function loads the baseline 2001 Land Cover raster, handles the 
#' Sinusoidal-to-WGS84 projection, and calculates the density (fraction) 
#' of each target class within a 0.15 degree (3x3 pixel) window.
#'
#' @param lc_file Path to the baseline Land Cover HDF file (e.g., 2001).
#' @param target_classes Vector of numeric IDs for the classes to process.
#' @param template_rast The 0.05 degree template raster for alignment.
#' @param study_geom The SpatVector defining West Africa.
#' @return A SpatRaster stack with one density layer per target class.
calculate_lc_density <- function(lc_file, target_classes, template_rast, study_geom) {
  library(terra)
  library(stringr)
  
  # --- 1. Load and Align the Baseline 2001 Raster ---
  LC_LAYER_NAME <- "LC_Type1"
  LC_FILL_VALUE <- 255
  
  file_df <- data.frame(path = lc_file) %>%
    mutate(Year = str_extract(basename(path), "(?<=A)\\d{4}")) %>%
    group_by(Year) %>%
    summarise(year_files = list(path), .groups = "drop") %>%
    arrange(Year)
  
  # --- 1. Load and Align Rasters ---
  lc_stack_list <- list()
  
  for (i in 1:nrow(file_df)) {
    current_year <- file_df$Year[i]
    tiles <- unlist(file_df$year_files[i])
    
    # A. Load tiles and select layer 1
    tile_rasters <- lapply(tiles, function(tile_path) {
      r <- terra::rast(tile_path)[LC_LAYER_NAME]
      return(r)
    })
    # Create a SpatRasterCollection
    sprc_tiles <- terra::sprc(tile_rasters)
    # Mosaic them (merge spatially)
    mosaic_year <- terra::mosaic(sprc_tiles, fun = "first") # 'first' is fine as overlaps match
    
    # B. Explicit CRS assignment
    if (is.na(crs(mosaic_year)) || crs(mosaic_year) == "") {
      # Standard MODIS Sinusoidal PROJ string
      crs(mosaic_year) <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    }
    
    # C. Crop to Study Area
    # Project vector to Sinusoidal
    geom_sinu <- terra::project(study_geom, crs(mosaic_year))
    lc_cropped <- terra::crop(mosaic_year, terra::ext(geom_sinu))
    
    # Project raster to WGS84 (template CRS)
    # method="near" for categorical data
    lc_wgs84 <- terra::project(lc_cropped, crs(template_rast), method = "near")
    
    # D. Resample to Template
    lc_aligned <- terra::resample(lc_wgs84, template_rast, method = "near")
    
    # E. Clean Fill Values (255 -> NA)
    lc_clean <- terra::subst(lc_aligned, LC_FILL_VALUE, NA)
    
    lc_stack_list[[i]] <- lc_clean
  }
  
  lc_base <- terra::rast(lc_stack_list)
  
  # --- 2. Calculate Focal Density for Target Classes ---
  
  density_list <- list()
  
  for (i in seq_along(target_classes)) {
    class_id <- target_classes[i]
    
    # A. Create Binary Map (1 = Class Present, 0 = Absent)
    # ifel is efficient for this
    binary_layer <- terra::ifel(lc_base == class_id, 1, 0)
    
    # B. Focal Calculation (3x3 window)
    # A 3x3 window on a 0.05 deg grid = 0.15 deg x 0.15 deg area.
    # Calculating the 'mean' of 1s and 0s gives the fraction (density).
    # w = 3 implies a 3x3 matrix of weights (all 1)
    density_layer <- terra::focal(binary_layer, w = 3, fun = "mean", na.rm = TRUE)
    
    # C. Name the layer
    # Naming convention: LC_[ID]_Density
    names(density_layer) <- paste0("LC_", class_id, "_Density")
    
    density_list[[i]] <- density_layer
  }
  
  # 3. Stack results
  return(terra::rast(density_list))
}
