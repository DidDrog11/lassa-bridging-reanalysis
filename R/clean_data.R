################################################################################
## SCRIPT: clean_data.R
##
## PURPOSE: Processes all raw data acquired in get_data.R. This script:
##          1. Assembles the finalized, aligned predictor stack (occ.covs).
##          2. Transforms multi-species rodent data into IMSOM-ready 3D arrays.
##
## INPUTS: All files saved to ./data/raw/ (rasters, Project_ArHa_database.rds, etc.)
## OUTPUTS: Saves the final predictor stack and the IMSOM input lists to ./data/processed/
################################################################################

# 1. Setup and Package Loading -------------------------------------------------

# Load core packages (defined in packages.R)
source(here::here("packages.R"))

# Load custom functions (defined in R folder)
source(here("R", "calculate_colwell_indices.R"))
source(here("R", "composite_ndvi.R"))
source(here("R", "calculate_lc_duration.R"))
source(here("R", "calculate_lc_density.R"))

# Define input/output directories
raw_dir <- here("data", "raw")
proc_dir <- here("data", "processed")
if (!dir.exists(proc_dir)) { dir.create(proc_dir, recursive = TRUE) }

# ---------------------------------------------------------------------------- #
# SECTION A.1: Define Study Grid Template (0.05 degree)
# ---------------------------------------------------------------------------- #
# This file.exists checks if we have already produced the final predictor stack
# It can be reproduced by removing the ! or running the code within the if manually
output_stack_path <- here(proc_dir, "final_predictor_stack.tif")
if(!file.exists(output_stack_path)) {
  # A.1.a. Read in raw GADM boundary files
  # GADM files are saved in the data/raw/gadm/ folder.
  gadm_files <- list.files(here(raw_dir, "gadm"), pattern = "\\_1_pk.rds$", full.names = TRUE)
  
  # Read all shapefiles into a list of SpatVectors
  list_of_spatvecs <- purrr::map(gadm_files, terra::vect)
  
  # Combine the individual country SpatVectors into a single object
  # This defines the final boundary of the study region (west_africa_ext) when ext() is used
  west_africa_ext <- do.call(rbind, list_of_spatvecs)
  
  # A.1.b. Load a high-resolution global raster to use for initial cropping
  raw_bio1 <- rast(here(raw_dir, "climate", "wc2.1_30s", "wc2.1_30s_bio_1.tif"))
  
  # Crop the high-resolution raster stack to the precise extent of the SpatVector boundary
  # This creates an intermediate raster that is used only to define the grid extent.
  cropped_reference <- crop(raw_bio1, west_africa_ext)
  
  # Define the final template SpatRaster at the required 0.05 degree resolution
  template_resolution <- 0.05
  
  # Create the definitive template raster object
  template_rast <- rast(extent = ext(cropped_reference), # Use the extent of the cropped data
                        resolution = template_resolution,
                        crs = "EPSG:4326")
  
  # Clean up large intermediate objects
  rm(raw_bio1, cropped_reference, list_of_spatvecs)
  
  # ---------------------------------------------------------------------------- #
  # SECTION A.2: Process Mean and Static Variables (Tmu, Pmu, Elev, Pop)
  # ---------------------------------------------------------------------------- #
  
  # A.2.a. Process WorldClim Proxies (Tmu, Pmu) ----------------------------------
  Tmu_proxy_wc <- rast(here(raw_dir, "climate", "wc2.1_30s", "wc2.1_30s_bio_1.tif")) |>
    crop(ext(west_africa_ext))
  Pmu_proxy_wc <- rast(here(raw_dir, "climate", "wc2.1_30s", "wc2.1_30s_bio_12.tif")) |>
    crop(ext(west_africa_ext))
  
  # 1. Unit Conversion (C, mm to mm/day)
  Tmu_final <- Tmu_proxy_wc
  Pmu_final <- Pmu_proxy_wc / 365.25 # Mean daily precip
  
  # 2. Resample and align to the study template
  # Use bilinear for continuous, smooth data like temperature and precipitation.
  Tmu_final <- resample(Tmu_final, template_rast, method = "bilinear")
  Pmu_final <- resample(Pmu_final, template_rast, method = "bilinear")
  
  names(Tmu_final) <- "Tmu"
  names(Pmu_final) <- "Pmu"
  
  # A.2.b. Process Elevation (Elev) ---------------------------------------------
  Elev_raw <- terra::rast(here(raw_dir, "elevation", "wc2.1_30s", "wc2.1_30s_elev.tif")) |>
    crop(ext(west_africa_ext))
  
  # Resample and align. Use bilinear interpolation for elevation.
  Elev_final <- resample(Elev_raw, template_rast, method = "bilinear")
  names(Elev_final) <- "Elev"
  
  # A.2.c. Process Population Density (Pop) -------------------------------------
  pop_dir <- here(raw_dir, "population_density") 
  pop_files <- lapply(list.files(pop_dir, pattern = "\\.tif$", full.names = TRUE), rast)
  pop_mosaic <- mosaic(sprc(pop_files))
  rm(pop_files)
  
  # Resample and align.
  # 'average' is safer for density aggregation to a coarser resolution.
  Pop_density <- pop_mosaic |> 
    terra::crop(west_africa_ext) |> 
    terra::resample(template_rast, method = "average")
  
  names(Pop_density) <- "Pop" # This goes into the model stack (Density)
  
  # Calculate Population count
  # Calculate area of each pixel in km2
  pixel_area <- terra::cellSize(Pop_density, unit = "km")
  # Count = Density * Area
  Pop_count <- Pop_density * pixel_area
  names(Pop_count) <- "Pop_Count_Raw"
  
  # Generate Urban Classification (Dijkstra Method)
  # Definition:
  # - City: Contiguous cells > 1,500 dens/km2 with Total Pop > 50,000
  # - Town: Contiguous cells > 300 dens/km2 with Total Pop > 5,000
  # - Rural: Everything else
  
  # Identify Candidate Pixels
  # Create masks for density thresholds
  dens_high <- terra::clamp(Pop_density, lower = 1500, values = FALSE) # NA if < 1500
  dens_mod  <- terra::clamp(Pop_density, lower = 300, values = FALSE)  # NA if < 300
  
  # Find Cities (High Density Clusters)
  # patches() groups connected pixels into unique IDs
  city_patches <- terra::patches(dens_high, directions = 8, zeroAsNA = TRUE)
  
  if (!all(is.na(values(city_patches)))) {
    # Sum population per patch
    city_pops <- terra::zonal(Pop_count, city_patches, fun = "sum")
    # Identify patch IDs that meet the 50k threshold
    valid_cities <- city_pops[city_pops[,2] > 50000, 1]
    # Create the City Raster (3)
    city_mask <- terra::match(city_patches, valid_cities, nomatch=NA)
    city_mask <- !is.na(city_mask) # TRUE where valid city
  } else {
    city_mask <- terra::rast(Pop_density)
    values(city_mask) <- 0
  }
  
  # Find Towns (Moderate Density Clusters)
  # Note: Towns include Cities spatially, so we classify "Towns and Cities" first
  town_patches <- terra::patches(dens_mod, directions = 8, zeroAsNA = TRUE)
  
  if (!all(is.na(values(town_patches)))) {
    town_pops <- terra::zonal(Pop_count, town_patches, fun="sum")
    valid_towns <- town_pops[town_pops[,2] > 5000, 1]
    town_mask <- terra::match(town_patches, valid_towns, nomatch=NA)
    town_mask <- !is.na(town_mask)
  } else {
    town_mask <- terra::rast(Pop_density)
    values(town_mask) <- 0
  }
  
  # Assemble Final Classification
  # Start with Rural (1)
  Urban_Class <- Pop_density
  names(Urban_Class) <- "Urban Class"
  values(Urban_Class) <- ifelse(is.na(values(Pop_density)), NA, 1)
  
  # Overlay Towns (2)
  Urban_Class[town_mask == 1] <- 2
  
  # Overlay Cities (3) - Cities supersede Towns
  Urban_Class[city_mask == 1] <- 3
  
  names(Urban_Class) <- "Urban_Class_Dijkstra"
  
  terra::writeRaster(Pop_count, here(proc_dir, "pop_count_2020_05deg.tif"), overwrite = TRUE)
  terra::writeRaster(Urban_Class, here(proc_dir, "pop_class_dijkstra_05deg.tif"), overwrite = TRUE)
  
  # A.2.d. Assemble Static/Mean Stack
  Static_Stack <- c(Tmu_final, Pmu_final, Elev_final, Pop_density)
  
  # ---------------------------------------------------------------------------- #
  # SECTION A.3: Precipitation and NDVI Layers
  # ---------------------------------------------------------------------------- #
  
  # A.3.a. Process CHIRPS ---------------------------------------------------
  chirps_dir <- here(raw_dir, "chirps")
  chirps_files <- list.files(chirps_dir, pattern = "\\.tif$", full.names = TRUE)
  
  # Total number of months in the stack (2001-01-01 to 2025-06-30)
  M_total <- length(chirps_files)
  
  P_stack_raw <- rast(chirps_files)
  
  # Set a threshold to capture all non-possible values.
  na_threshold <- -100
  # Define the rules matrix for substitution:
  # Column 1: Lower bound (inclusive)
  # Column 2: Upper bound (inclusive)
  # Column 3: Replacement value (NA = Not Available)
  rules_matrix <- matrix(c(-Inf,      # Lower bound: Negative infinity
                           na_threshold, # Upper bound: The threshold (-100)
                           NA         # Replacement value: NA
  ), ncol = 3)
  P_stack_clean <- classify(P_stack_raw, rcl = rules_matrix, include.lowest = TRUE)
  
  # Assuming the layers are ordered chronologically.
  # This should be true if the list.files() result is alphabetically ordered by filename (YYYYMM).
  
  # I. Mean (Pmu) and Variability (Pcv)
  # We calculate Pmu using the full 2001-2025 data
  Pmu_recalc <- app(P_stack_clean, fun = "mean", na.rm = TRUE) / 30 # Convert to monthly
  Pcv <- app(P_stack_clean, fun = function(x) { 
    # Coefficient of Variation (CV) = (Standard Deviation / Mean)
    sd_val <- sd(x, na.rm = TRUE)
    mean_val <- mean(x, na.rm = TRUE)
    return(sd_val / mean_val)
  })
  names(Pcv) <- "Pcv"
  
  # II. Extremes (Pmin, Pmax)
  Pmin <- app(P_stack_clean, fun = "min", na.rm = TRUE, names = "Pmin")
  Pmax <- app(P_stack_clean, fun = "max", na.rm = TRUE, names = "Pmax")
  
  # III. Duration of Low Precipitation (Pdur)
  # Action: Calculate the number of months per year below 1 mm/day (or 30.4 mm/month)
  # The original paper used 1 mm/day. Since these are monthly sums, we use a monthly threshold.
  # Threshold: 1 mm/day * ~30 days/month = ~30 mm/month
  P_low_flag <- app(P_stack_clean, fun = function(x) { x < 30 })
  # Sum the low months over the entire time series, then divide by the number of years.
  P_low_flag_sum <- app(P_low_flag, fun = "sum", na.rm = TRUE)
  # Calculate the average number of low-flag months per year (M_total / 12)
  Pdur <- P_low_flag_sum / (M_total / 12)
  names(Pdur) <- "Pdur"
  
  # IV. Colwell's Indices (Pc, Pm) - Contingency and Constancy
  P_colwell_stack <- app(P_stack_clean, fun = calculate_colwell_indices)
  names(P_colwell_stack) <- c("Pc", "Pm")
  # Extract the individual layers from the resulting stack
  Pc <- P_colwell_stack[["Pc"]]
  Pm <- P_colwell_stack[["Pm"]]
  
  # Combine all precipitation layers
  P_indices_stack <- c(Pmu_recalc, Pcv, Pmin, Pmax, Pdur, Pc, Pm)
  target_names <- c("Pmu", "Pcv", "Pmin", "Pmax", "Pdur", "Pc", "Pm")
  # Apply the names to the SpatRaster object
  names(P_indices_stack) <- target_names
  
  # A.3.b. NDVI -------------------------------------------------------------
  ndvi_dir <- here(raw_dir, "ndvi")
  ndvi_files <- list.files(ndvi_dir, pattern = "\\.hdf$", full.names = TRUE)
  N_stack_clean_path <- here(proc_dir, "ndvi_clean_monthly_stack.tif")
  
  # Action: Spatially aggregate MOD and MYD products into a single monthly composite.
  # We assume the layers are named such that they can be temporally ordered and grouped 
  if (file.exists(N_stack_clean_path)) {
    N_stack_clean <- terra::rast(N_stack_clean_path)
  } else {
    # 1. Spatially aggregate MOD and MYD products into a single monthly composite.
    N_stack_clean <- composite_ndvi(file_list = ndvi_files, 
                                    template_rast = template_rast, 
                                    study_area_geom = west_africa_ext)
    
    # 2. Save the result for future runs
    terra::writeRaster(N_stack_clean, N_stack_clean_path, overwrite = TRUE)
  }
  
  # The raw composite still has values like 8000. We multiply by 0.0001.
  N_stack_scaled <- N_stack_clean * 0.0001
  
  Y_total <- M_total / 12
  
  # I. Mean (Nmu) and Variability (Ncv)
  Nmu <- terra::app(N_stack_scaled, fun = "mean", na.rm = TRUE)
  names(Nmu) <- "Nmu"
  
  Ncv <- terra::app(N_stack_scaled, fun = function(x) { 
    mean_val <- mean(x, na.rm = TRUE)
    # 1. Check for NA
    if (is.na(mean_val)) { 
      return(NA) 
    }
    # 2. Check for Negative or Zero Mean
    # CV is only ecologically meaningful when the mean is positive.
    if (mean_val <= 0) { 
      return(NA) # Set to NA to avoid negative CV and division by zero
    }
    # Calculate CV if the mean is valid and positive
    sd_val <- sd(x, na.rm = TRUE)
    return(sd_val / mean_val)
  })
  names(Ncv) <- "Ncv"
  
  # II. Extremes (Nmin, Nmax)
  Nmin <- terra::app(N_stack_scaled, fun = "min", na.rm = TRUE, names = "Nmin")
  Nmax <- terra::app(N_stack_scaled, fun = "max", na.rm = TRUE, names = "Nmax")
  
  # III. Duration of Low NDVI (Ndur)
  # Threshold: NDVI < 0.5. Normalizing by actual total years.
  N_low_flag <- terra::app(N_stack_scaled, fun = function(x) { x < 0.5 })
  Ndur <- terra::app(N_low_flag, fun = "sum", na.rm = TRUE) / Y_total
  names(Ndur) <- "Ndur"
  
  # IV. Colwell's Indices (Nc, Nm) - Contingency and Constancy
  # Using the custom function defined in R/calculate_colwell_indices.R
  N_colwell_stack <- terra::app(N_stack_scaled, fun = calculate_colwell_indices)
  
  # Extract and rename the layers
  names(N_colwell_stack) <- c("Nc", "Nm")
  Nc <- N_colwell_stack[["Nc"]]
  Nm <- N_colwell_stack[["Nm"]]
  
  
  # 3. Combine all NDVI layers
  N_indices_stack <- c(Nmu, Ncv, Nmin, Nmax, Ndur, Nc, Nm)
  target_names <- c("Nmu", "Ncv", "Nmin", "Nmax", "Ndur", "Nc", "Nm")
  # Apply the names to the SpatRaster object
  names(N_indices_stack) <- target_names

  # A.3.c. Nighttime Lights (NTL) -------------------------------------------
  ntl_dir <- here(raw_dir, "ntl")
  ntl_files <- list.files(ntl_dir, pattern = "median_masked.*\\.tif$", full.names = TRUE)
  
  NTL_raw <- terra::rast(ntl_files[1])
  # 1. Crop/Resample to template
  # NTL is ~500m resolution, so we aggregate (average) to 0.05 deg (~5km)
  # Average preserves the "intensity density" of the pixel
  NTL_aligned <- NTL_raw |> 
    terra::crop(west_africa_ext) |> 
    terra::resample(template_rast, method = "average")
    
  # 2. Log-Transformation
  # NTL data is extremely right-skewed (Cities are bright, rural is zero).
  # We apply log1p(x) = log(x + 1) to linearise the relationship for models.
  NTL_log <- log1p(NTL_aligned)
    
  names(NTL_log) <- "NTL"
  
  Static_Stack_NTL <- c(Static_Stack, NTL_log)
  
  # ---------------------------------------------------------------------------- #
  # SECTION A.4: Land Cover Stability and Density Layers
  # ---------------------------------------------------------------------------- #
  # A.4.b. Stability Analysis (Find stable classes) -------------------------
  lc_dir <- here(raw_dir, "landcover")
  lc_files <- list.files(lc_dir, pattern = "\\.hdf$", full.names = TRUE)
  
  # Define path for the cached stable classes table
  stable_classes_path <- here(proc_dir, "stable_lc_classes.csv")
  
  if (file.exists(stable_classes_path)) {
    stable_classes_df <- read.csv(stable_classes_path)
    
  } else {
    # Run the duration calculation function
    # This will handle the Sinusoidal projection, cropping, and transition matrix
    stable_classes_df <- calculate_lc_duration(
      lc_files = lc_files, 
      study_geom = west_africa_ext, 
      template_rast = template_rast
    )
    
    # Save the list of stable classes for future runs
    write.csv(stable_classes_df, stable_classes_path, row.names = FALSE)
  }
  
  # Save the list of stable classes for reference
  write.csv(stable_classes_df, here(proc_dir, "stable_lc_classes.csv"), row.names = FALSE)
  message(paste("Identified", nrow(stable_classes_df), "stable Land Cover classes."))
  
  # A.4.b Perform Density Analysis ------------------------------------------
  lc_2001_file <- lc_files[grep("A2001", basename(lc_files))]
  
  # Extract the stable class IDs
  target_class_ids <- stable_classes_df$LC_ID
  
  # 2. FORCE INCLUSION of Ecologically Critical Classes
  # LC 12 = Croplands (Critical for Mastomys)
  # LC 2  = Evergreen Forest (Critical for exclusion/Malacomys)
  # LC 14 = Crop Mosaic (Useful if pure crop is rare)
  forced_classes <- c(12, 2, 14)
  
  final_target_ids <- unique(c(target_class_ids, forced_classes))
  final_target_ids <- sort(final_target_ids)
  
  # Run the density calculation
  LC_density_stack <- calculate_lc_density(
    lc_file = lc_2001_file,
    target_classes = final_target_ids,
    template_rast = template_rast,
    study_geom = west_africa_ext
  )
  

  # A.5: Final Predictor Stack Assembly -----------------------------------
  # Combine all components: 
  # 1. Static/Mean (Tmu, Pmu, Elev, Pop, NTL)
  # 2. Precipitation Indices (Pcv, Pmin, Pmax, Pdur, Pc, Pm, etc.)
  # 3. NDVI Indices (Nmu, Ncv, Nmin, Nmax, Ndur, Nc, Nm)
  # 4. Land Cover Densities
  
  Static_Stack_aligned <- terra::resample(Static_Stack_NTL, template_rast, method = "near") |>
    select(-"Pmu")
  P_indices_aligned <- terra::resample(P_indices_stack, template_rast, method = "near")
  N_indices_aligned <- terra::resample(N_indices_stack, template_rast, method = "near")
  LC_density_aligned <- terra::resample(LC_density_stack, template_rast, method = "near")
  
  # Combine all component stacks
  Final_Predictor_Stack <- c(Static_Stack_aligned, 
                             P_indices_aligned, 
                             N_indices_aligned, 
                             LC_density_aligned)
  
  # Save the final stack
  # Note: we later log transform some layers before scaling the entire stack
  output_stack_path <- here(proc_dir, "final_predictor_stack.tif")
  terra::writeRaster(Final_Predictor_Stack, output_stack_path, overwrite = TRUE)
} else {
  Final_Predictor_Stack <- rast(output_stack_path)
}

# ---------------------------------------------------------------------------- #
# SECTION B: Rodent Data Structuring (IMSOM Inputs)
# ---------------------------------------------------------------------------- #

# B.1. Load and Standardise Data Sources ----------------------------------
# 1. Load Project ArHa Database
arha_path <- here(raw_dir, "Project_ArHa_database.rds")
if (!file.exists(arha_path)) stop("Project ArHa RDS not found.")
arha_db <- readRDS(arha_path)

# 2. Define Target Scope
# Countries (ISO3)
wa_poly_list <- lapply(list.files(path = here("data", "raw", "gadm"), pattern = "_1_pk", full.names = T), function(x)
  vect(x))
west_africa_ext <- do.call(rbind, wa_poly_list)
target_countries <- unique(west_africa_ext$GID_0)

# Species (Strict JSDM Set)
# This is a subset of the species_to_obtain defined in get_data.R
target_species <- tibble(rodent = c(
  "Mastomys natalensis",      # Primary LASV Reservoir (Clade A-I target)
  "Mastomys erythroleucus",   # Congeneric competitor & potential secondary reservoir
  "Rattus rattus",            # Invasive dominant competitor (Urban/Peridomestic exclusion)
  "Mus musculus",             # Invasive domestic competitor
  "Praomys rostratus",        # Common native species (often co-occurs in peridomestic zones)
  "Lophuromys sikapusi",      # Native species with historical LASV antibody detections
  "Arvicanthis niloticus",    # Common savanna species (habitat overlap)
  "Malacomys edwardsi"        # Forest specialist (Ecological control/contrast)
))

target_species <- left_join(target_species, arha_db$host |>
                              distinct(host_species, gbif_id),
                            by = c("rodent" = "host_species"))

if(!file.exists(here(proc_dir, "rodent_data_mapped.rds"))) {
  # 3. Process ArHa Data
  # A. Initial Filter and Date Parsing
  rodent_arha_base <- arha_db$host |>
    filter(iso3c %in% target_countries) |>
    filter(gbif_id %in% target_species$gbif_id) |>
    filter(!is.na(latitude) & !is.na(longitude)) |>
    mutate(source = "ArHa_PA",
           date_obj = lubridate::ymd(start_date),
           year = lubridate::year(date_obj),
           replicate_id = format(date_obj, "%Y-%m"), # Group by Month-Year
           observed = ifelse(number_of_hosts > 0, 1, 0))
  
  # B. Calculate Effort Statistics
  # Calculate mean/sd only from the data where effort is explicitly recorded.
  valid_effort_idx <- !is.na(rodent_arha_base$trap_nights_clean)
  valid_effort_values <- rodent_arha_base$trap_nights_clean[valid_effort_idx]
  
  # Log-transform to handle the range (58 to 12,000)
  log_effort <- log(valid_effort_values + 1)
  mean_log_effort <- mean(log_effort)
  sd_log_effort <- sd(log_effort)
  
  # C. Finalise ArHa structure
  rodent_arha <- rodent_arha_base |>
    mutate(effort_raw = trap_nights_clean,
           effort_known = ifelse(!is.na(effort_raw), 1, 0),
           # Standardized Effort Variable
           # If known: (Log(x) - Mean) / SD
           # If unknown: Set to 0 (The "Mean" value). This neutralises the slope.
           effort_scaled = ifelse(!is.na(effort_raw), (log(effort_raw + 1) - mean_log_effort) / sd_log_effort, 0)) |>
    dplyr::select(record_id = host_record_id, study_id, species = host_species,
                  lat = latitude, lon = longitude, year, replicate_id,
                  observed, effort_scaled, effort_known, source)
  
  # D. Additional urban rodent data
  urban_lit_path <- here(raw_dir, "rodent_urban_literature.csv")
  rodent_urban_clean <- tibble()
  
  if(file.exists(urban_lit_path)) {
    # Load and clean
    urban_raw <- read.csv(urban_lit_path, stringsAsFactors = FALSE)
    
    rodent_urban_processed <- urban_raw |>
      filter(species %in% target_species$rodent) |>
      mutate(effort_numeric = as.numeric(ifelse(trap_nights == "not_reported", NA, trap_nights)),
             # Determine Source Type
             # If effort is known, treat as High Quality (Source 1). If unknown, Opportunistic (Source 2).
             source_type = ifelse(!is.na(effort_numeric), "UrbanLit_HQ", "UrbanLit_PO"),
             # Replicate ID: We default to Month 01
             replicate_id = paste0(year, "-01"), 
             # Observation Status
             # For HQ data: number_of_hosts = n_detected
             number_of_hosts = n_detected,
             observed = ifelse(n_detected > 0, 1, 0),
             # IDs
             record_id = paste0("Urb_", row_number()),
             study_id = paste0("Lit_", source) # Using the 'source' column (DOI/Ref)
      ) |>
      dplyr::select(record_id, study_id, species, lat, lon, year, replicate_id,
                    number_of_hosts, observed, effort_numeric, source = source_type)
  }
  
  urban_hq <- rodent_urban_processed |>
    filter(source == "UrbanLit_HQ") |>
    select(-observed) |>
    mutate(trap_nights_clean = effort_numeric) # Match ArHa column name
  
  # E. West African Rodents dataset
  wa_rds_path <- here(raw_dir, "WA_rodents_database.rds")
  rodent_wa_hq <- tibble()
  
  if(exists("wa_db")) {
    
    # Helper: Parse DMS strings (e.g., "09_50_27")
    parse_dms <- function(x) {
      # If NA or empty, return NA
      if(is.na(x) || x == "") return(NA_real_)
      # Split by underscore
      parts <- as.numeric(unlist(strsplit(x, "_")))
      if(length(parts) != 3) return(NA_real_)
      # D + M/60 + S/3600
      return(parts[1] + parts[2]/60 + parts[3]/3600)
    }
    
    parse_utm <- function(utm_col, zone_default=28) {
      lats <- numeric(length(utm_col)) * NA
      lons <- numeric(length(utm_col)) * NA
      # Identify non-NA UTMs
      valid_idxs <- which(!is.na(utm_col) & utm_col != "")
      if(length(valid_idxs) > 0) {
        for(i in valid_idxs) {
          parts <- unlist(strsplit(utm_col[i], "_"))
          if(length(parts) == 3) {
            zone_str <- gsub("[A-Z]", "", parts[1])
            zone <- as.numeric(zone_str)
            easting <- as.numeric(parts[2])
            northing <- as.numeric(parts[3])
            epsg_code <- 32600 + zone
            pt <- sf::st_point(c(easting, northing))
            sfc <- sf::st_sfc(pt, crs = epsg_code)
            pt_wgs <- sf::st_transform(sfc, 4326)
            coords <- sf::st_coordinates(pt_wgs)
            lons[i] <- coords[1]
            lats[i] <- coords[2]
          }
        }
      }
      return(list(lat = lats, lon = lons))
    }
    
    utm_results <- parse_utm(wa_db$UTM_coordinates)
    wa_db$lat_utm <- utm_results$lat
    wa_db$lon_utm <- utm_results$lon
    
    # Cleaning Pipeline
    rodent_wa_clean <- wa_db |>
      # 1. Clean Taxonomy
      # Remove records not identified to species (contains "-" or "sp." or "spp.")
      filter(!grepl("-|sp\\.|spp\\.", species)) |>
      mutate(
        # Construct scientific name: "Genus species"
        scientificName = paste(str_to_title(genus), str_to_lower(species))
      ) |>
      # Filter to target species
      filter(scientificName %in% target_species$rodent) |>
      # 2. Clean Dates (EventDate)
      mutate(
        # Year: First 4 digits of year_trapping
        year_clean = as.numeric(str_extract(year_trapping, "^\\d{4}")),
        # Fallback Year: Extract from unique_id (fc_2010...) and subtract 2 (publication delay)
        year_fallback = as.numeric(str_extract(unique_id, "\\d{4}")) - 2,
        year_final = coalesce(year_clean, year_fallback),
        # Month: Extract first 3 chars, map to number
        # "Apr-Sep" -> "Apr" -> 04
        month_str = str_to_title(str_sub(month_trapping, 1, 3)),
        month_num = match(month_str, month.abb), # Built-in R constant month.abb
        # Construct Replicate ID (YYYY-MM)
        # Drop records with no valid month (spring/autumn/NA)
        month_final = ifelse(is.na(month_num), "01", sprintf("%02d", month_num)),
        replicate_id = paste0(year_final, "-", month_final)
      ) |>
      filter(!is.na(year_final)) |>
      # 3. Clean Coordinates
      rowwise() |>
      rowwise() |>
      mutate(
        lat_dms = parse_dms(latitude_DMS_N),
        lon_dms = parse_dms(longitude_DMS_W) * -1, 
        lat_final = coalesce(as.numeric(latitude_D_N), lat_utm, lat_dms),
        lon_final = coalesce(as.numeric(longitude_D_E), lon_utm, lon_dms)
      ) |>
      ungroup() |>
      filter(!is.na(lat_final) & !is.na(lon_final)) |>
      # 4. Clean Metrics (Counts and Effort)
      mutate(
        # Number detected
        observed = ifelse(as.numeric(number) > 0, 1, 0),
        effort_val = suppressWarnings(as.numeric(trap_nights)),
        source = "WA_Rodents_PA",
        record_id = paste0("WA_", row_number()),
        study_id = unique_id
      ) |>
      dplyr::select(record_id, study_id, species = scientificName, 
                    lat = lat_final, lon = lon_final, year = year_final, 
                    replicate_id, observed, effort_val, source)
    
    # Assign to the output object
    rodent_wa_hq <- rodent_wa_clean
    message(paste("WA Rodents Cleaned:", nrow(rodent_wa_hq), "records retained."))
  }
  
  # --- 4. DEDUPLICATION ---
  # Goal: Create a unified Source 1 (HQ) dataset without duplicates.
  # Trust Order: ArHa > UrbanLit > WA_Rodents
  
  # Function to generate spatial-temporal keys
  get_key <- function(df) {
    # Round coords to 2dp (~1.1km) to catch slight precision differences
    paste(round(df$lat, 2), round(df$lon, 2), df$year, df$species, sep = "_")
  }
  
  # A. Use ArHa as the core
  source1_core <- rodent_arha
  core_keys <- get_key(source1_core)
  
  # B. Filter WA_Rodents and urban rodents against Core
  if(nrow(rodent_wa_hq) > 0) {
    wa_keys <- get_key(rodent_wa_hq)
    urban_keys <- get_key(urban_hq)
    # Keep only unique records
    rodent_wa_clean <- rodent_wa_hq[!wa_keys %in% core_keys, ]
    rodent_urban_clean <- urban_hq[!urban_keys %in% core_keys, ]
    message(paste("WA Rodents: Loaded", nrow(rodent_wa_hq), "| Retained", nrow(rodent_wa_clean), "after deduplication.\n"),
            paste("Urban Rodents: Loaded", nrow(urban_hq), "| Retained", nrow(rodent_urban_clean), "after deduplication."))
  } else {
    rodent_wa_clean <- tibble()
  }
  
  rodent_hq_raw <- bind_rows(rodent_arha, rodent_urban_clean, rodent_wa_clean)
  
  # Consolidate and Rescale Effort
  rodent_hq_clean <- rodent_hq_raw |>
    mutate(
      effort_raw_combined = coalesce(trap_nights_clean, effort_numeric, effort_val)
    )
  
  # Calculate Global Scaling Statistics (Log-Scale)
  # Only use records where effort is known
  valid_effort <- na.omit(rodent_hq_clean$effort_raw_combined)
  mean_log_effort <- mean(log(valid_effort + 1))
  sd_log_effort <- sd(log(valid_effort + 1))
  
  message(paste("Global Effort Mean (Log):", round(mean_log_effort, 2)))
  message(paste("Global Effort SD (Log):", round(sd_log_effort, 2)))
  
  # Final Cleanup
  source1_master <- rodent_hq_clean |>
    mutate(
      effort_scaled = ifelse(!is.na(effort_raw_combined), 
                             (log(effort_raw_combined + 1) - mean_log_effort) / sd_log_effort, 
                             0), # 0 = Mean effort (imputed for missing)
      effort_known = ifelse(!is.na(effort_raw_combined), 1, 0)
    )  |>
    dplyr::select(record_id, study_id, species, lat, lon, year, replicate_id, 
                  observed, effort_scaled, effort_known, source)
  
  # 4. Process & Filter Opportunistic Data (Original + GBIF + Urban rodents)
  # Helper to format opportunistic sources uniformly
  format_opportunistic <- function(df, src_name) {
    df |>
      mutate(source = src_name, observed = 1, effort_scaled = 0, effort_known = 0) |>
      dplyr::select(record_id, study_id, species, lat, lon, year, replicate_id, 
                    observed, effort_scaled, effort_known, source)
  }
  
  # A. Load Original Paper Data
  orig_path <- here(raw_dir, "mastomys_presences_original.csv")
  rodent_orig <- tibble()
  if(file.exists(orig_path)) {
    rodent_orig <- read.csv(orig_path) |>
      filter(Species == "Mastomys natalensis") |>
      mutate(record_id = paste0("Orig_", row_number())) |>
      group_by(Reference) |>
      mutate(study_id = paste0("Orig_Study_", cur_group_id())) |>
      ungroup() |>
      mutate(year = as.numeric(stringr::str_extract(Year, "\\d{4}")),
             replicate_id = paste0(year, "-01")) |>
      rename(species = Species, lat = Latitude, lon = Longitude) |>
      format_opportunistic("Original_PO")
  }
  
  # B. Load GBIF Data
  gbif_path <- here(raw_dir, "gbif_data.csv")
  rodent_gbif <- tibble()
  
  if (file.exists(gbif_path)) {
    # 1. Load Raw Data
    # read_tsv is safer for GBIF data than read.csv
    gbif_raw <- read.csv(gbif_path, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE) 
    
    # 2. Stricter species IDs
    # Species that are difficult to ID morphologically or are commensal
    difficult_id_species <- c("Mus musculus", 
                              "Mastomys natalensis", 
                              "Mastomys erythroleucus", 
                              "Rattus rattus")
    
    # Accepted evidence types for difficult species
    high_quality_basis <- c("PRESERVED_SPECIMEN", "MATERIAL_SAMPLE", "HUMAN_OBSERVATION")
    
    # Institutions with suspected data issues (e.g. uniform counts, unlikely ecology)
    suspicious_institutions <- c("Laboratoire des Sciences Forestières (LSF/UAC)")
    
    # 3. Filter and Process
    rodent_gbif <- gbif_raw |>
      # A. Taxonomy and Geography Filters
      filter(speciesKey %in% target_species$gbif_id) %>% # Match GBIF IDs
      filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) |>
      filter(year >= 2001) |>
      # B. QC Filter: Remove Suspicious Institutions
      filter(!institutionCode %in% suspicious_institutions) |>
      # C. QC Filter: Strict Evidence for Difficult Species
      # Logic: If species is "Easier" (e.g. Arvicanthis), keep all records.
      filter(case_when(species %in% difficult_id_species ~ basisOfRecord %in% high_quality_basis,
                       TRUE ~ TRUE)) |>
      # D. Formatting
      mutate(record_id = as.character(gbifID),
             study_id = paste0("GBIF_", datasetKey),
             replicate_id = paste0(year, "-01"), # Default Replicate
             species = species, 
             lat = decimalLatitude, 
             lon = decimalLongitude) |>
      format_opportunistic("GBIF_PO")
  } else {
    warning("GBIF data file not found.")
  }
  
  urban_po <- rodent_urban_processed |>
    filter(source == "UrbanLit_PO") |>
    filter(observed == 1) |> # Keep only presences for PO stream
    mutate(effort_scaled = 0, effort_known = 0) |>
    select(-number_of_hosts)
  
  # C. Combine Opportunistic Data
  rodent_opp <- bind_rows(rodent_orig, rodent_gbif, urban_po)
  
  # 4. Generate pseudo-absences for opportunistic data
  # We need background points (0s) to contrast with the GBIF/Original presences (1s).
  # These represent the "Available Environment".
  
  # Target number of background points (2000, but less than the recommended 10,000 due to computational requirements)
  TARGET_BACKGROUND_N <- 2000
  
  # Sample random points from the template raster (ensures valid land pixels)
  # as.points=TRUE returns a SpatVector
  pseudo_points <- terra::spatSample(Final_Predictor_Stack, 
                                     size = TARGET_BACKGROUND_N, 
                                     method = "random", 
                                     na.rm = TRUE, 
                                     as.points = TRUE)
  
  # Convert to Dataframe and format to match rodent_opp_clean
  pseudo_df <- data.frame(terra::geom(pseudo_points)) |>
    dplyr::select(x, y) |>
    mutate(# Metadata to match the Opportunistic structure
      record_id = paste0("Pseudo_", row_number()),
      study_id = "Background_Sampling",
      species = "Pseudo_Absence", # Placeholder species
      lat = y, 
      lon = x,
      year = 2015, # Midpoint year (or random sampling 2001-2025)
      replicate_id = "Background_Rep",
      observed = 0, # Observation Data
      # Covariates (Opportunistic = 0 effort, 0 quality)
      effort_scaled = 0,
      effort_known = 0,
      source = "Opportunistic_PO") |>
    dplyr::select(-x, -y)
  
  # Combine Presences (rodent_opp_clean) with Pseudo-Absences (pseudo_df)
  rodent_opp_final <- bind_rows(rodent_opp, pseudo_df)
  
  keys_s1 <- get_key(source1_master)
  keys_s2 <- get_key(rodent_opp_final)
  
  # Identify Source 2 records that duplicate a Source 1 record
  dupe_indices <- which(keys_s2 %in% keys_s1)
  
  if (length(dupe_indices) > 0) {
    source2_clean <- rodent_opp_final[-dupe_indices, ]
    message(paste("Final Safety Check: Removed", length(dupe_indices), 
                  "records from Source 2 that duplicated High-Quality Source 1 data."))
  } else {
    source2_clean <- rodent_opp_final
    message("Final Safety Check: No cross-source duplicates found.")
  }
  
  # Use the cleand Source 2
  all_rodent_data <- bind_rows(source1_master, source2_clean) %>%
    mutate(species = as.factor(species))
  
  # Check final counts
  table(all_rodent_data$source)
  
  message(paste("Source 2 Final Count:", nrow(rodent_opp_final), 
                "(", nrow(rodent_opp), "Presences +", nrow(pseudo_df), "Background)"))
  
  # B.2. Spatial Aggregation (Map to 0.05 deg Grid) -------------------------
  if (!exists("Final_Predictor_Stack")) {
    Final_Predictor_Stack <- terra::rast(here(proc_dir, "final_predictor_stack.tif"))
  }
  
  # Create SpatVector
  rodent_vect <- terra::vect(all_rodent_data, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # Extract Cell Numbers
  cells_extracted <- terra::extract(Final_Predictor_Stack, rodent_vect, cells = TRUE, ID = FALSE)
  
  # Attach Cell ID
  all_rodent_data$cell_id <- cells_extracted$cell
  
  # Filter invalid points
  all_rodent_data_clean <- all_rodent_data |>
    filter(!is.na(cell_id)) |>
    mutate(species = as.character(species))
  
  # Save Mapped Data
  write_rds(all_rodent_data_clean, here(proc_dir, "rodent_data_mapped.rds"))
} else {
  
  all_rodent_data_clean <- read_rds(here(proc_dir, "rodent_data_mapped.rds"))
  
}

# B.3. Generate IMSOM Inputs (3D Arrays & Index Lists) --------------------
# --- 1. Create Master Site List and Extract Covariates (occ.covs) ---
# The model needs a dataframe of predictors for every unique cell_id.

# Identify all unique sites (J)
master_cell_ids <- sort(unique(all_rodent_data_clean$cell_id))

if(!file.exists(here(proc_dir, "final_predictor_stack_scaled.tif"))) {
  # --- Selective Log-Transformation ---
  # Load/Scale the Predictor Stack if not already done
  if (!exists("Final_Predictor_Stack")) {
    Final_Predictor_Stack <- terra::rast(here(proc_dir, "final_predictor_stack.tif"))
  }
  
  # A. Transformations (Log-Transform Skewed Vars)
  # Only run if not already transformed in memory
  if (terra::minmax(Final_Predictor_Stack["Pop"])[2] > 20) {
    skewed_vars <- c("Pop", "LC_13_Density")
    for (var in skewed_vars) {
      if (var %in% names(Final_Predictor_Stack)) {
        Final_Predictor_Stack[[var]] <- terra::app(Final_Predictor_Stack[[var]], fun = log1p)
      }
    }
  }
  # Global Scaling (Z-scores)
  Final_Predictor_Stack_Scaled <- terra::scale(Final_Predictor_Stack, center = TRUE, scale = TRUE)
  scaled_stack_path <- here(proc_dir, "final_predictor_stack_scaled.tif")
  terra::writeRaster(Final_Predictor_Stack_Scaled, scaled_stack_path, overwrite = TRUE)
} else {
  Final_Predictor_Stack_Scaled <- rast(here(proc_dir, "final_predictor_stack_scaled.tif"))
}

# Extract predictors from the stack for these specific cells
# This ensures row 1 of occ_covs matches site index 1
occ_covs_scaled <- terra::extract(Final_Predictor_Stack_Scaled, master_cell_ids)

# Add the cell_id as a column for linking
occ_covs_scaled$cell_id <- master_cell_ids

# Filter for NAs (Sites with no environmental data, e.g. ocean edge)
valid_sites_idx <- complete.cases(occ_covs_scaled)
if (sum(!valid_sites_idx) > 0) {
  occ_covs_scaled <- occ_covs_scaled[valid_sites_idx, ]
  master_cell_ids <- occ_covs_scaled$cell_id
  # Filter rodent data to match valid sites
  all_rodent_data_clean <- all_rodent_data_clean %>% 
    filter(cell_id %in% master_cell_ids)
}

# 2. Define Dimensions and Helper Function
species_levels <- sort(unique(all_rodent_data_clean$species))
species_levels <- species_levels[species_levels != "Pseudo_Absence"]
N_species <- length(species_levels)

# Consolidate Visits (Many rows -> One row per Visit)
# This function collapses the long data so we have one row per Site/Replicate.
# It aggregates the observed species into a list-column.
consolidate_visits <- function(df) {
  df |>
    group_by(cell_id, replicate_id) |>
    summarise(
      effort_val = first(effort_scaled),
      quality_val = first(effort_known),
      visit_source = first(source),
      # Observation: List of species found (observed == 1) at this visit
      found_species = list(unique(species[observed == 1])),
      .groups = "drop"
    )
}

# 3. Apply consolidation to sources
# Source 1: ArHa (High Quality)
source1_raw <- all_rodent_data_clean |> filter(grepl("PA|HQ", source))
source1_consolidated <- consolidate_visits(source1_raw)

# Source 2: Opportunistic (GBIF + Original + Pseudo)
# Regex matches "GBIF", "Original", or "Opportunistic" (Pseudo-absences)
source2_raw <- all_rodent_data_clean |> filter(grepl("PO", source))
source2_consolidated <- consolidate_visits(source2_raw)

message(paste("Source 1 Visits:", nrow(source1_consolidated)))
message(paste("Source 2 Visits:", nrow(source2_consolidated)))

# Master Site List (J): All unique grid cells involved
J_total <- length(master_cell_ids)

message(paste("Global Dimensions: J =", J_total, "| N =", N_species))

# Build Source 1 Array (ArHa - High Quality)
# A. Identify Sites and Dimensions for Source 1
s1_site_ids <- sort(unique(source1_consolidated$cell_id))
J_s1 <- length(s1_site_ids)

# Calculate Max Replicates (K) for Source 1
reps_s1 <- source1_consolidated |> count(cell_id) |> pull(n)
K_s1 <- max(reps_s1)
message(paste("Source 1 Dimensions: Sites =", J_s1, "| Max Reps =", K_s1))

# B. Initialize Arrays with Defaults
# y: Default is NA (representing "Visit did not happen")
y_s1 <- array(NA, dim = c(N_species, J_s1, K_s1))
# Covariates: Default is 0 (Neutral mean for standardized covariates)
# This ensures padding replicates have valid numeric values for the model matrix
effort_s1  <- array(0, dim = c(J_s1, K_s1))
quality_s1 <- array(0, dim = c(J_s1, K_s1))

# C. Fill Arrays Loop
for (j in 1:J_s1) {
  # Get the specific cell_id
  current_site <- s1_site_ids[j]
  
  # Get visits for this site (chronological)
  site_visits <- source1_consolidated |> 
    filter(cell_id == current_site) |> 
    arrange(replicate_id)
  
  # Loop through actual visits
  for (k in 1:nrow(site_visits)) {
    visit <- site_visits[k, ]
    # 1. Fill Covariates
    effort_s1[j, k]  <- visit$effort_val
    quality_s1[j, k] <- visit$quality_val
    # 2. Fill Observations (y)
    # Since visit happened, initialize all species to 0 (Not Detected)
    y_s1[, j, k] <- 0
    # 3. Overwrite Presences
    # 'found_species' is a list-column of names
    present_spp <- unlist(visit$found_species)
    
    if (length(present_spp) > 0) {
      spp_indices <- match(present_spp, species_levels)
      spp_indices <- spp_indices[!is.na(spp_indices)]
      if (length(spp_indices) > 0) {
        y_s1[spp_indices, j, k] <- 1
      }
    }
  }
}
# D. Create Site Map (Index into Master List)
sites_s1_map <- match(s1_site_ids, master_cell_ids)

# Build Source 2 Array (Opportunistic: Original + GBIF + Pseudo)
# A. Identify Sites and Dimensions
s2_site_ids <- sort(unique(source2_consolidated$cell_id))
J_s2 <- length(s2_site_ids)

reps_s2 <- source2_consolidated |> count(cell_id) |> pull(n)
K_s2 <- max(reps_s2) # Should be 1, but code allows flexibility
message(paste("Source 2 Dimensions: Sites =", J_s2, "| Max Reps =", K_s2))

# B. Initialize Arrays
y_s2 <- array(NA, dim = c(N_species, J_s2, K_s2))
# Covariates are 0 (Intercept only model)
effort_s2  <- array(0, dim = c(J_s2, K_s2))
quality_s2 <- array(0, dim = c(J_s2, K_s2))

# C. Fill Arrays Loop
for (j in 1:J_s2) {
  current_site <- s2_site_ids[j]
  
  site_visits <- source2_consolidated |> 
    filter(cell_id == current_site) |> 
    arrange(replicate_id)
  
  for (k in 1:nrow(site_visits)) {
    visit <- site_visits[k, ]
    
    # 1. Covariates (Stay 0)
    
    # 2. Fill Observations (y)
    # Initialize all to 0 (Not Detected / Background)
    y_s2[, j, k] <- 0
    
    # 3. Overwrite Presences
    present_spp <- unlist(visit$found_species)
    if (length(present_spp) > 0) {
      spp_indices <- match(present_spp, species_levels)
      spp_indices <- spp_indices[!is.na(spp_indices)]
      
      if (length(spp_indices) > 0) {
        y_s2[spp_indices, j, k] <- 1
      }
    }
  }
}
# D. Create Site Map
sites_s2_map <- match(s2_site_ids, master_cell_ids)

# Compile Final List for spOccupancy
imsom_input_list <- list(
  # Observations
  y = list(y_s1, y_s2),
  # Site Covariates
  occ.covs = occ_covs_scaled,
  # Detection Covariates (List of Lists)
  det.covs = list(list(effort = effort_s1, quality = quality_s1), # Source 1
                  list(effort = effort_s2, quality = quality_s2)  # Source 2 (Dummy vars)
  ),
  # Site Indices
  sites = list(sites_s1_map, sites_s2_map),
  # Species List
  species_names = species_levels,
  # Coordinates
  coords = terra::xyFromCell(Final_Predictor_Stack, master_cell_ids)
)

# 5. Save
saveRDS(imsom_input_list, here(proc_dir, "imsom_input_list.rds"))

# ---------------------------------------------------------------------------- #
# SECTION C: Epidemiological Data Structuring
# ---------------------------------------------------------------------------- #

# C.1. Rodent Pathogen Data (ArHa Primary) --------------------------------
lasv_rodent_arha <- arha_db$pathogen |>
  left_join(arha_db$host, by = "host_record_id") |>
  filter(host_species == "Mastomys natalensis",
         pathogen_species_cleaned == "Mammarenavirus lassaense",
         !is.na(latitude), !is.na(longitude),
         number_tested > 0) |>
  # Aggregate to Site Level
  group_by(latitude, longitude) |>
  summarise(n_tested = sum(number_tested, na.rm = TRUE),
            n_pos = sum(number_positive, na.rm = TRUE),
            .groups = "drop") |>
  mutate(outcome = case_when(n_pos > 0 ~ 1, # Define Outcome: 1 = Found, 0 = Not Found (if effort sufficient)
                             n_tested >= 1 ~ 0,
                             TRUE ~ NA_real_),
         source = "ArHa_pathogen") |> # Ambiguous
  filter(!is.na(outcome))

# 2. Process New Sierra Leone Data (Line List to Site Aggregation)
# Load file (Send this to Git at somepoint)
sl_raw_path <- "C:/Users/ucbtds4/R_Repositories/data_for_pharos/data/rodent_eastern_sl_2023-08-11.csv"
if(file.exists(sl_raw_path)) {
  sl_data <- readr::read_csv(sl_raw_path, show_col_types = FALSE) |>
    filter(`Host species` == "Mastomys natalensis",
           grepl("Lassa", `Detection target`, ignore.case = TRUE),
           !is.na(Latitude), !is.na(Longitude)) |>
    mutate(is_positive = ifelse(`Detection outcome` == "positive", 1, 0),
           # Extract Village Code (e.g., "LAL" from "1_LAL_001")
           village_code = stringr::str_split_i(`Sample ID`, "_", 2)) |>
    # Aggregate to Village Level
    group_by(village_code) |>
    summarise(latitude = mean(Latitude, na.rm = TRUE), # Centroid
              longitude = mean(Longitude, na.rm = TRUE),
              n_tested = n(),
              n_pos = sum(is_positive, na.rm = TRUE),
              source = "SL_pathogen",
              .groups = "drop")
} else {
  sl_data <- tibble()
}

# 3. Process Basinski Rodent Data
basinski_path <- here(raw_dir, "human_seroprevalence_basinski.csv")
if(file.exists(basinski_path)) {
  basinski_rodent_raw <- readr::read_csv(basinski_path) |>
    filter(Species == "natalensis",
           !is.na(Latitude), !is.na(Longitude)) |>
    select(latitude = Latitude, longitude = Longitude, 
           n_pos = NumPosAb, n_tested = NumTestAb) |>
    mutate(source = "Basinski_pathogen") |>
    # Handle potential NAs in their counts
    filter(!is.na(n_tested) & n_tested > 0)
  
  # --- SPATIAL DUPLICATE CHECK ---
  # We must remove Basinski points that are within 500m of ArHa or New SL points
  # to avoid double counting the same historical studies.
  
  # Combine current master set
  current_master <- bind_rows(lasv_rodent_arha, sl_data)
  
  if(nrow(current_master) > 0 && nrow(basinski_rodent_raw) > 0) {
    v_master <- terra::vect(current_master, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    v_basinski <- terra::vect(basinski_rodent_raw, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    
    # Find neighbours within (500m)
    nearby <- terra::nearby(v_basinski, v_master, distance = 500)
    
    # Get indices of Basinski points that have a neighbour in Master
    duplicates_idx <- unique(nearby[,1])
    
    if(length(duplicates_idx) > 0) {
      basinski_rodent_clean <- basinski_rodent_raw[-duplicates_idx, ]
    } else {
      basinski_rodent_clean <- basinski_rodent_raw
    }
  } else {
    basinski_rodent_clean <- basinski_rodent_raw
  }
  
} else {
  basinski_rodent_clean <- tibble()
}

# 4. Combine and Define Outcome
# Outcome Logic:
# 1 = Found (at least 1 positive)
# 0 = Not Found (if at least 1 tested)
# NA = Ambiguous (Not found but low effort)
lasv_rodent <- bind_rows(lasv_rodent_arha, sl_data, basinski_rodent_clean) |>
  mutate(outcome = case_when(n_pos > 0 ~ 1,
                             n_tested >= 1 ~ 0,
                             TRUE ~ NA_real_ )) |>
  filter(!is.na(outcome)) |>
  select(-village_code)

# C.2. Human Seroprevalence Data (Calibration) ----------------------------
# A. Basinski Data (Baseline)
human_base <- readr::read_csv(here(raw_dir, "human_seroprevalence_basinski.csv")) |>
  filter(Species == "sapiens") |>
  mutate(year_explicit = as.numeric(stringr::str_extract(Year, "\\d{4}")),
         pub_year = as.numeric(stringr::str_extract(Source, "\\d{4}")),
         year_final = case_when(!is.na(year_explicit) ~ year_explicit,
                           !is.na(pub_year) ~ pub_year - 2,
                           TRUE ~ NA_real_)) |>
  dplyr::select(lat = Latitude, lon = Longitude, 
                n_test = NumTestAb, n_pos = NumPosAb, 
                ref = Source, year = year_final) |>
  filter(!is.na(year)) |>
  mutate(source = "Basinski_pathogen")

# B. New Data (Additional Sources)
new_human_path <- here(raw_dir, "human_seroprevalence_updated.csv")
if(file.exists(new_human_path)) {
  human_new <- readr::read_delim(new_human_path, delim = ";") |>
    filter(Human_Random_Survey == TRUE) |>
    select(lat = Latitude, lon = Longitude, n_test = NumTestAb, n_pos = NumPosAb, ref = Source, year = Year) |>
    mutate(source = "Updated_pathogen",
           lat = as.numeric(lat),
           lon = as.numeric(lon))
  
  # Bind them together
  all_human_data <- bind_rows(human_base, human_new)
} else {
  all_human_data <- human_base
  warning("Additional human seroprevalence file not found. Using Basinski baseline only.")
}

# B.2. Synthetic Urban Anchors
# Rationale: Add high-confidence absence points in hyper-urban centres to 
# constrain the model against extrapolating risk into high-NTL zones.
urban_anchors <- data.frame(id = c("Syn_Lagos_Island", "Syn_Abidjan_Plateau", "Syn_Accra_Osu", "Syn_Conakry_Dixinn", "Syn_Freetown_CBD"),
                            lat = c(6.4522, 5.33973, 5.55191, 9.6035, 8.461),
                            lon = c(3.40406, -4.02512, -0.20944, -13.62901, -13.229),
                            n_pos = c(0, 0, 0, 0, 0),
                            n_test = c(500, 500, 500, 500, 500),
                            ref = "Constraint_Urban_Core",
                            year = 2025,
                            source = "Synthetic_Constraint")

# C. Clean Human Data
# Ensure valid coordinates and positive N
all_human_data <- bind_rows(all_human_data, urban_anchors) |>
  filter(!is.na(lat) & !is.na(lon) & n_test > 0)

# C.3. Spatial Blocking (Anti-Circularity w/ retention for sparse  --------
# Ensure we have the Human Data to block against
if (exists("all_human_data") && nrow(all_human_data) > 0) {
  # 1. Convert to vect objects
  human_vect <- vect(all_human_data, geom = c("lon", "lat"), crs = "EPSG:4326")
  rodent_vect <- vect(lasv_rodent, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  
  # 2. Find overlaps using nearby()
  # Returns a matrix with 2 columns: [id_rodent, id_human]
  # distance is in meters (5000m = 5km)
  nearby_matches <- terra::nearby(rodent_vect, human_vect, distance = 5000)
  
  # 3. Define Data-Sparse Countries (Where we CANNOT afford to lose rodent data)
  # Ghana, Benin, Togo, Cote d'Ivoire are sparse for Lassa testing.
  # Nigeria/Sierra Leone/Guinea are relatively data-rich.
  sparse_countries <- rbind(human_vect, rodent_vect) |>
    terra::intersect(west_africa_ext |>
                           select(COUNTRY)) |>
    as_tibble() |>
    group_by(COUNTRY) |>
    summarise(n = n()) |>
    filter(n <= 20) |>
    pull(COUNTRY)
  
  # 4. Apply Rescue/Block Logic
  # The country of each rodent point involved in a match
  rodent_countries <- terra::intersect(rodent_vect, west_africa_ext[, "COUNTRY"])
  
  rodent_drop_indices <- c()
  human_drop_indices <- c()
  
  if (nrow(nearby_matches) > 0) {
    for (i in 1:nrow(nearby_matches)) {
      r_idx <- nearby_matches[i, 1]
      h_idx <- nearby_matches[i, 2]
      
      # Get Country of the Rodent point
      # Note: Handling potential NA if point is slightly offshore
      ctry <- rodent_countries$COUNTRY[r_idx] 
      if (ctry %in% sparse_countries) {
        # RESCUE: Data is scarce here. 
        # Keep Rodent (for Training), Drop Human (from Validation)
        human_drop_indices <- c(human_drop_indices, h_idx)
      } else {
        # STANDARD: Data is rich. 
        # Drop Rodent (prevent circularity), Keep Human (prioritize validation)
        rodent_drop_indices <- c(rodent_drop_indices, r_idx)
      }
    }
  }
  
  # 5. Filter Datasets
  
  # A. Filter Rodents
  if (length(rodent_drop_indices) > 0) {
    rodent_drop_indices <- unique(rodent_drop_indices)
    lasv_rodent_clean <- lasv_rodent[-rodent_drop_indices, ]
    message(paste("Blocked", length(rodent_drop_indices), "rodent sites in data-rich regions."))
  } else {
    lasv_rodent_clean <- lasv_rodent
  }
  
  # B. Filter Humans (Update the validation set)
  if (length(human_drop_indices) > 0) {
    human_drop_indices <- unique(human_drop_indices)
    all_human_data_final <- all_human_data[-human_drop_indices, ]
    message(paste("Removed", length(human_drop_indices), "human sites in sparse regions to allow rodent training."))
  } else {
    all_human_data_final <- all_human_data
  }
  
} else {
  warning("Human data not found. Skipping spatial blocking.")
  lasv_rodent_clean <- lasv_rodent
  all_human_data_final <- all_human_data
}

message(paste("Final Rodent Training N:", nrow(lasv_rodent_clean)))
message(paste("Final Human Validation N:", nrow(all_human_data_final)))


# C.4. Prepare and Save Training Data -------------------------------------
# Extract Predictors for Rodent Data (for BRT)
stack_path <- here(proc_dir, "final_predictor_stack_scaled.tif")

# Garbage collection to free up RAM before the heavy extract operation
gc() 

if (file.exists(stack_path)) {
  pred_stack_clean <- terra::rast(stack_path)
} else {
  stop("Predictor stack not found. Section C.4 failed.")
}

train_vect <- terra::vect(lasv_rodent_clean, geom = c("longitude", "latitude"), crs = "EPSG:4326")
train_covs <- terra::extract(pred_stack_clean, train_vect, ID = FALSE)

# Combine predictors with outcome
pathogen_training_data <- cbind(lasv_rodent_clean, train_covs) |>
  na.omit() 

# Save Outputs
saveRDS(pathogen_training_data, here(proc_dir, "pathogen_training_data.rds"))
write.csv(all_human_data, here(proc_dir, "human_calibration_data.csv"), row.names = FALSE)

# D.1 Load and Standardise Human Case Data --------------------------------

moore_raw_path <- here(raw_dir, "lassa_cases_moore.csv")

if (file.exists(moore_raw_path)) {
  # Load raw
  moore_data <- readr::read_csv(moore_raw_path)
  
  # Filter for Validation Window (2012-2022)
  # We focus on Admin 1 and 2
  val_window_start <- 2012
  val_window_end   <- 2022
  n_years_window <- val_window_end - val_window_start + 1
  
  moore_clean <- moore_data |>
    filter(ADM0_name %in% target_countries,
           YEAR >= val_window_start & YEAR <= val_window_end,
           Lowest_gadm_level_entered %in% c(1, 2)) |>
    mutate(Admin_Level = Lowest_gadm_level_entered,
           Location_Code = ifelse(Admin_Level == 2, ADM2_code, ADM1_code),
           # Fix encoding if possible, or just rely on Code
           Location_Name = ifelse(Admin_Level == 2, ADM2_name, ADM1_name)) |> 
    group_by(ADM0_name, Admin_Level, Location_Code, Location_Name) |>
    summarise(Total_Confirmed = sum(CASES, na.rm = TRUE),
              Data_Start_Year = min(YEAR, na.rm = TRUE),
              Data_End_Year = max(YEAR, na.rm = TRUE),
              .groups = "drop") |>
    mutate(Annual_Cases = Total_Confirmed / n_years_window)
  
  # Save Admin 2 Data (District/LGA)
  write.csv(filter(moore_clean, Admin_Level == 2), 
            here(proc_dir, "validation_cases_admin2.csv"), row.names = FALSE)
  
  # Save Admin 1 Data (State/Province)
  write.csv(filter(moore_clean, Admin_Level == 1), 
            here(proc_dir, "validation_cases_admin1.csv"), row.names = FALSE)
}



# E. Output focus sites ---------------------------------------------------
# Define in clean_data.R or load in every results script
standardised_cities <- data.frame(name = c(
  # Tier 1: Large cities / National Capitals (High NTL, Concrete Core)
  "Lagos", "Abidjan", "Accra", "Conakry", "Freetown",
  # Tier 2: Regional Cities (Medium NTL, Peri-urban dominant)
  "Bouake", "Tamale", "Ogbomosho", 
  # Tier 3: Endemic Towns
  "Kenema", "Ekpoma", "N'Zerekore", "Jos"),
  type = c(rep("Large City", 5),
           rep("Regional City", 3),
           rep("Endemic Town", 4)),
  # Coordinates centred on the "Urban Core"
  lat = c(6.4549, 5.3204, 5.5560, 9.5350 + (1 * 0.05), 8.4844,   # Large cities
          7.6903, 9.4075, 8.1333,                   # Regional
          7.8767, 6.7432, 7.7562, 9.92              # Towns
          ),
  lon = c(3.4246, -4.0161, -0.1969, -13.6800 + (3 * 0.05), -13.2344, # Large cities
          -5.0357, -0.8534, 4.2667,                     # Regional
          -11.1875, 6.1390, -8.8179, 8.89               # Towns
          ))

write_rds(standardised_cities, here(proc_dir, "standardised_cities.rds"))
