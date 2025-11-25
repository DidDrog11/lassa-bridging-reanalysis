# R/composite_ndvi.R

#' Creates a single monthly NDVI composite stack from raw MODIS HDF files
#'
#' Uses the {raster} package as a bridge to force data type stabilization 
#' before processing in {terra}.
#'
#' @param file_list A vector of HDF file paths (MOD/MYD).
#' @param template_rast The 0.05 degree reference SpatRaster for final alignment.
#' @param study_area_geom The SpatVector defining the West Africa study area boundaries.
#' @return A SpatRaster object containing one layer for each unique month.
composite_ndvi <- function(file_list, template_rast, study_area_geom) {
  # Ensure necessary packages are loaded in the calling script
  library(terra)
  library(stringr)
  library(dplyr)
  library(raster)
  
  # --- Define Constants ---
  NDVI_LAYER <- "CMG 0.05 Deg Monthly NDVI\""
  NDVI_SCALE_FACTOR <- 0.0001
  RAW_VALID_MIN <- -2000 # Valid raw integer range min
  RAW_VALID_MAX <- 10000 # Valid raw integer range max
  NDVI_FILL_VALUE <- -3000 # Fill value from documentation
  
  # 1. Prepare File Mapping
  file_df <- data.frame(path = file_list) |>
    mutate(date_str = str_extract(basename(path), "A[0-9]{7}"),
           year = str_sub(date_str, 2, 5),
           doy = as.numeric(str_sub(date_str, 6, 8)),
           date = as.Date(paste0(year, "-", doy), format = "%Y-%j"),
           year_month = format(date, "%Y%m")) |>
    group_by(year_month) %>%
    summarise(monthly_files = list(path), .groups = "drop") %>%
    arrange(year_month)
  
  # 2. Process and Composite Each Month
  composite_list <- list()
  study_extent <- terra::ext(study_area_geom)
  
  for (i in 1:nrow(file_df)) {
    month_files <- unlist(file_df$monthly_files[i])
    individual_layers <- list()
    
    for (j in 1:length(month_files)) {
      file <- month_files[j]
      
      # --- 1. Read and Select Layer using Raster Bridge ---
      
      # A. Get the full DSN path string (Subdataset 1)
      sds_metadata <- terra::describe(file)
      dsn_line <- sds_metadata[str_detect(sds_metadata, "SUBDATASET_1_NAME=")]
      # Extract everything after the equals sign
      full_dsn_path <- stringr::str_extract(dsn_line, "(?<=SUBDATASET_1_NAME=).*") 
      
      # B. Load using raster::raster() to force numeric type coercion
      rast_layer_bridge <- raster::raster(full_dsn_path, 
                                          # Tell 'raster' the specific scale factor
                                          scale = 0.0001, 
                                          # Tell 'raster' the specific Fill Value to convert to NA during read
                                          na.value = -3000)
      
      # C. Convert to a SpatRaster
      raw_ndvi_layer <- terra::rast(rast_layer_bridge) 
      
      # Define the CRS (WGS 84 to match GADM/template)
      terra::crs(raw_ndvi_layer) <- "EPSG:4326"
      
      # --- 2. Processing: Crop, Clean, and Scale ---
      # 1. CROP: Crop the raw layer to the study extent
      ndvi_cropped <- terra::crop(raw_ndvi_layer, study_area_geom)
      
      # 2. CLEAN: Remove erroneous fill values and outliers
      raw_values <- terra::values(ndvi_cropped)
      numeric_values <- suppressWarnings(as.numeric(raw_values))
      numeric_values[is.nan(numeric_values)] <- NA
      numeric_values[numeric_values == NDVI_FILL_VALUE] <- NA
      # Scale remaining values
      numeric_values = numeric_values * NDVI_SCALE_FACTOR
      # Remove non valid values
      numeric_values[numeric_values < RAW_VALID_MIN | numeric_values > RAW_VALID_MAX] <- NA
      
      # Step 2b: Reapply the values
      terra::values(ndvi_cropped) <- numeric_values
      
      # 4. ALIGN: Resample to the final 0.05 deg template
      individual_layers[[j]] <- terra::resample(ndvi_cropped, template_rast, method = "bilinear")
    }
    
    # 3. COMPOSITE: Combine MOD and MYD layers for the current month (using median)
    if (length(individual_layers) > 0) {
      monthly_stack_resampled <- terra::rast(individual_layers)
      composite_layer <- terra::app(monthly_stack_resampled, fun = "median", na.rm = TRUE)
      names(composite_layer) <- file_df$year_month[i]
      composite_list[[i]] <- composite_layer
    }
  }
  
  # 4. Combine all monthly composites into one stack
  return(terra::rast(composite_list))
}