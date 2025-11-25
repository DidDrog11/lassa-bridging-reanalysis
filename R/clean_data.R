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

# A.1.a. Read in raw GADM boundary files
# GADM files are saved in the data/raw/gadm/ folder.
gadm_files <- list.files(here(raw_dir, "gadm"), pattern = "\\.rds$", full.names = TRUE)

# Read all shapefiles into a list of SpatVectors
list_of_spatvecs <- purrr::map(gadm_files, terra::vect)

# Combine the individual country SpatVectors into a single object
# This defines the final boundary of the study region (west_africa_ext)
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
Pop_final <- pop_mosaic |> 
  crop(west_africa_ext) |>
  resample(template_rast, method = "average")

names(Pop_final) <- "Pop"


# A.2.d. Assemble Static/Mean Stack
Static_Stack <- c(Tmu_final, Pmu_final, Elev_final, Pop_final)

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

# ---------------------------------------------------------------------------- #
# SECTION A.4: Land Cover Stability and Density Layers
# ---------------------------------------------------------------------------- #

# --- A.4.a. Stability Analysis (Find stable classes) ---
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


# --- A.4.b. Focal Density Calculation ---
# Extract the stable class IDs
target_class_ids <- stable_classes_df$LC_ID

# Run the density calculation
LC_density_stack <- calculate_lc_density(
  lc_file = lc_2001_file,
  target_classes = target_class_ids,
  template_rast = template_rast,
  study_geom = west_africa_ext
)

# ---------------------------------------------------------------------------- #
# SECTION A.5: Final Predictor Stack Assembly
# ---------------------------------------------------------------------------- #
# Combine all components: 
# 1. Static/Mean (Tmu, Pmu, Elev, Pop)
# 2. Precipitation Indices (Pcv, Pmin, Pmax, Pdur, Pc, Pm, etc.)
# 3. NDVI Indices (Nmu, Ncv, Nmin, Nmax, Ndur, Nc, Nm)
# 4. Land Cover Densities (6 layers)

Static_Stack_aligned <- terra::resample(Static_Stack, template_rast, method = "near") |>
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
output_stack_path <- here(proc_dir, "final_predictor_stack.tif")
terra::writeRaster(Final_Predictor_Stack, output_stack_path, overwrite = TRUE)


