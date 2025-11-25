# R/calculate_lc_duration.R

#' Calculates the average duration (T_i) for each Land Cover class based on yearly transitions.
#'
#' Uses pure {terra} to load 25 years of MODIS Land Cover (MCD12Q1), 
#' tracks pixel-by-pixel transitions, and calculates the stability of each class.
#'
#' @param lc_files A vector of file paths for the 25 yearly MCD12Q1 rasters (2001-2025).
#' @param study_geom The SpatVector defining the West Africa study area boundaries.
#' @param template_rast The 0.05 degree template raster for alignment.
#' @return A data frame listing the stable classes (T_i >= 20 years) and their IDs/Durations.
calculate_lc_duration <- function(lc_files, study_geom, template_rast) {
  library(terra)
  library(dplyr)
  library(stringr)
  
  # --- Constants ---
  LC_LAYER_NAME <- "LC_Type1"
  LC_FILL_VALUE <- 255
  
  file_df <- data.frame(path = lc_files) %>%
    # Extract Year from filename (e.g., MCD12Q1.A2001001...)
    # Pattern: A followed by 4 digits (Year)
    mutate(Year = str_extract(basename(path), "(?<=A)\\d{4}")) %>%
    group_by(Year) %>%
    summarise(year_files = list(path), .groups = "drop") %>%
    arrange(Year)
  
  # --- 1. Load and Align Yearly Rasters ---
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
  
  # Combine into a 25-layer stack
  lc_stack <- terra::rast(lc_stack_list)
  names(lc_stack) <- paste0("Year_", 2001:(2000 + terra::nlyr(lc_stack)))
  
  # --- 2. Transition Tally ---
  
  # Extract pixels to a matrix (Pixels x Years)
  lc_matrix <- terra::values(lc_stack, dataframe = FALSE)
  
  # Identify all unique Land Cover IDs present in the region (ignoring NA/Fill)
  # 255 is the fill value for MODIS Land Cover
  unique_lc_ids <- sort(unique(as.vector(lc_matrix)))
  unique_lc_ids <- unique_lc_ids[!is.na(unique_lc_ids) & unique_lc_ids != 255]
  
  # Initialize Tally Data Frame
  tally <- data.frame(
    LC_ID = unique_lc_ids,
    Appeared_Count = 0,
    Stayed_Count = 0
  )
  
  # Iterate through transitions (Year t -> Year t+1)
  n_years <- ncol(lc_matrix)
  
  for (t in 1:(n_years - 1)) {
    curr_col <- lc_matrix[, t]
    next_col <- lc_matrix[, t + 1]
    
    # Loop through each class found in the region
    for (i in 1:nrow(tally)) {
      id <- tally$LC_ID[i]
      
      # 1. Where did the class appear in year t?
      appeared_mask <- (curr_col == id)
      
      # 2. Where did it stay the same in year t+1?
      # (Must calculate mask on valid pixels only)
      stayed_mask <- appeared_mask & (next_col == id)
      
      # Sum ignoring NAs (NA counts as 0)
      tally$Appeared_Count[i] <- tally$Appeared_Count[i] + sum(appeared_mask, na.rm = TRUE)
      tally$Stayed_Count[i] <- tally$Stayed_Count[i] + sum(stayed_mask, na.rm = TRUE)
    }
  }
  
  
  # --- 3. Calculate Metrics (M_ii and T_i) ---
  results <- tally %>%
    filter(Appeared_Count > 0) %>%
    mutate(
      # Probability of remaining in class i (M_ii)
      M_ii = Stayed_Count / Appeared_Count,
      
      # Probability of changing (C_i)
      C_i = 1 - M_ii,
      
      # Expected Duration in Years (T_i = 1 / C_i)
      # Handle perfect stability (C_i = 0) by assigning a high max duration (e.g. 100 years)
      Duration_Years = ifelse(C_i <= 0.001, 100, 1 / C_i)
    )
  
  
  # --- 4. Filter and Label ---
  # Standard IGBP Legend
  lc_names <- c(
    "0" = "Water_Bodies",
    "1" = "Evergreen_Needleleaf_Forest",
    "2" = "Evergreen_Broadleaf_Forest",
    "3" = "Deciduous_Needleleaf_Forest",
    "4" = "Deciduous_Broadleaf_Forest",
    "5" = "Mixed_Forest",
    "6" = "Closed_Shrublands",
    "7" = "Open_Shrublands",
    "8" = "Woody_Savannas",
    "9" = "Savannas",
    "10" = "Grasslands",
    "11" = "Permanent_Wetlands",
    "12" = "Croplands",
    "13" = "Urban_and_Built-up",
    "14" = "Cropland_Natural_Vegetation_Mosaic",
    "15" = "Snow_and_Ice",
    "16" = "Barren_or_Sparsely_Vegetated",
    "17" = "Water_Bodies"
  )
  
  results$LC_Name <- lc_names[as.character(results$LC_ID)]
  
  # Filter for stability >= 20 years
  stable_classes <- results %>%
    filter(Duration_Years >= 20) %>%
    select(LC_ID, LC_Name, Duration_Years, M_ii)
  
  return(stable_classes)
}