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
# Note: we later log transform some layers before scaling the entire stack
output_stack_path <- here(proc_dir, "final_predictor_stack.tif")
terra::writeRaster(Final_Predictor_Stack, output_stack_path, overwrite = TRUE)

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

# 3. Process ArHa Data
# A. Initial Filter and Date Parsing
rodent_arha_base <- arha_db$host |>
  filter(iso3c %in% target_countries) |>
  filter(gbif_id %in% target_species$gbif_id) |>
  filter(!is.na(latitude) & !is.na(longitude)) |>
  mutate(source_type = "ArHa_PA",
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
                observed, effort_scaled, effort_known, source = source_type)

# 4. Process & Filter Opportunistic Data (Original + GBIF)
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
# --- B. Load GBIF Data (Source: GBIF_PO) ---
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
  high_quality_basis <- c("PRESERVED_SPECIMEN", "MATERIAL_SAMPLE")
  
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
    # Logic: If species is "Difficult", MUST be Specimen/Sample. 
    #        If species is "Easier" (e.g. Arvicanthis), keep all records.
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

# C. Combine Opportunistic Data
rodent_opp <- bind_rows(rodent_orig, rodent_gbif)

# D. Spatial Duplicate Filtering
if (nrow(rodent_opp) > 0) {
  # Convert both datasets to SpatVectors
  v_arha <- terra::vect(rodent_arha, geom = c("lon", "lat"), crs = "EPSG:4326")
  v_opp  <- terra::vect(rodent_opp, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # Find nearby points: returns a matrix of indices [id_opp, id_arha]
  # distance = 5000 meters (5 km)
  nearby_matrix <- terra::nearby(v_opp, v_arha, distance = 5000)
  
  # We only care about matches where the year is also the same
  # Extract the matching pairs
  opp_indices <- nearby_matrix[,1]
  arha_indices <- nearby_matrix[,2]
  
  # Compare years for the spatial matches
  match_year <- rodent_opp$year[opp_indices] == rodent_arha$year[arha_indices]
  
  # Identify Opportunistic IDs that are spatio-temporal duplicates
  ids_to_remove <- unique(opp_indices[match_year])
  
  # Filter the opportunistic dataset
  # We keep rows that are NOT in the removal list
  if (length(ids_to_remove) > 0) {
    rodent_opp_clean <- rodent_opp[-ids_to_remove, ]
  } else {
    rodent_opp_clean <- rodent_opp
  }
} else {
  rodent_opp_clean <- tibble()
}

# 5. Combine All Clean Data
all_rodent_data <- bind_rows(rodent_arha, rodent_opp_clean) %>%
  mutate(species = as.factor(species))

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
  filter(!is.na(cell_id))

# Save Mapped Data
write_rds(all_rodent_data_clean, here(proc_dir, "rodent_data_mapped.rds"))

# B.3. Generate IMSOM Inputs (3D Arrays & Index Lists) --------------------
# --- 1. Create Master Site List and Extract Covariates (occ.covs) ---
# The model needs a dataframe of predictors for every unique cell_id.

# Identify all unique sites (J)
master_cell_ids <- sort(unique(all_rodent_data_clean$cell_id))

# --- Selective Log-Transformation ---
# Skewed Variables
skewed_vars <- c("Pop", "LC_13_Density")

# Apply Log1p (Log(x + 1)) Before Scaling
for (var in skewed_vars) {
  if (var %in% names(Final_Predictor_Stack)) {
    # Apply function to the specific layer
    Final_Predictor_Stack[[var]] <- terra::app(Final_Predictor_Stack[[var]], fun = log1p)
  }
}

# Global Scaling (Z-scores)
Final_Predictor_Stack_Scaled <- terra::scale(Final_Predictor_Stack, center = TRUE, scale = TRUE)

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
# ----------------------------------------------------------------
species_levels <- sort(unique(all_rodent_data_clean$species))
species_levels <- species_levels[species_levels != "Pseudo_Absence"]
N_species <- length(species_levels)

# Function to create 3D array for a specific data source
format_source_array <- function(df, source_label, all_sites, all_spp) {
  # Filter data for this source
  # Use Regex to catch "ArHa" or "Opportunistic"
  df_sub <- df %>% filter(grepl(source_label, source))
  if (nrow(df_sub) == 0) return(NULL)
  # Identify sites present in this source
  # (Subset of the master list)
  source_site_ids <- sort(unique(df_sub$cell_id))
  J_source <- length(source_site_ids)
  # Determine K (Max Replicates)
  # Count visits per site
  visits <- df_sub %>% 
    group_by(cell_id) %>% 
    summarise(n_reps = n_distinct(replicate_id))
  K_max <- max(visits$n_reps)
  # Initialize Arrays
  # Y: N x J x K
  y_arr <- array(NA, dim = c(length(all_spp), J_source, K_max))
  # Det Covs: J x K (Effort)
  effort_arr <- array(NA, dim = c(J_source, K_max))
  quality_arr <- array(NA, dim = c(J_source, K_max))
  # Fill Arrays
  for (j in 1:J_source) {
    site <- source_site_ids[j]
    # Get records for this site
    site_recs <- df_sub[df_sub$cell_id == site, ]
    # Identify unique replicates (time points)
    reps <- sort(unique(site_recs$replicate_id))
    for (k in 1:length(reps)) {
      rep_id <- reps[k]
      # Get records for this specific visit
      visit_data <- site_recs[site_recs$replicate_id == rep_id, ]
      # 1. Fill Detection Covariates (Site x Rep)
      # Take the first value (effort should be constant for the visit)
      effort_arr[j, k] <- visit_data$effort_scaled[1]
      quality_arr[j, k] <- visit_data$effort_known[1]
      # 2. Fill Detection History (Species x Site x Rep)
      # Default: 0 (Not Detected) for ALL modeled species
      y_arr[, j, k] <- 0
      # Find which species were actually observed (1)
      present_spp <- visit_data |>
        filter(observed == 1) %>% pull(species)
      # Set 1s for present species
      if (length(present_spp) > 0) {
        spp_idx <- match(present_spp, all_spp)
        # Ignore NAs (e.g. if Pseudo_Absence was passed)
        spp_idx <- spp_idx[!is.na(spp_idx)]
        y_arr[spp_idx, j, k] <- 1
      }
    }
  }
  # Create Index Map (Where does this source fit in the Master Site List?)
  site_map <- match(source_site_ids, all_sites)
  return(list(y = y_arr,
              det_covs = list(effort = effort_arr, quality = quality_arr),
              sites = site_map))
}


# 3. Generate Inputs for Both Sources
# ----------------------------------------------------------------

# Source 1: ArHa (High Quality, Multi-Replicate)
source1 <- format_source_array(all_rodent_data_clean, "ArHa", master_cell_ids, species_levels)
# Source 2: Opportunistic (Presence-Only, Single Replicate)
source2 <- format_source_array(all_rodent_data_clean, "Original_PO|GBIF_PO", master_cell_ids, species_levels)

# 4. Compile Final List for spOccupancy
# ----------------------------------------------------------------
imsom_input_list <- list(
  # Observations
  y = list(source1$y, source2$y),
  # Site Covariates (Master DataFrame)
  occ.covs = occ_covs_scaled,
  # Detection Covariates (List of Lists)
  det.covs = list(source1$det_covs, source2$det_covs),
  # Site Indices (List of Vectors)
  sites = list(source1$sites, source2$sites),
  # Species List (Used for labelling output later)
  species_names = species_levels,
  # Coordinates (For spatial random effects later)
  coords = terra::xyFromCell(Final_Predictor_Stack, master_cell_ids)
)

# 5. Final checks
# Check for any non-finite coordinates
any_bad_coords <- any(!is.finite(imsom_input_list$coords))
if (any_bad_coords) {
  stop("CRITICAL: 'coords' contains NA, NaN, or Inf values. Check xyFromCell.")
} else {
  message("Coordinates are clean.")
}

# Species sparsity check
# Sum observations across sites and replicates for Source 1
obs_s1 <- apply(imsom_input_list$y[[1]], 1, sum, na.rm = TRUE)
# Sum observations across sites and replicates for Source 2
obs_s2 <- apply(imsom_input_list$y[[2]], 1, sum, na.rm = TRUE)

# Combine
total_obs <- obs_s1 + obs_s2
names(total_obs) <- imsom_input_list$species_names

print(total_obs)

if (any(total_obs < 5)) {
  warning("Some species have < 5 detections total. Model convergence is unlikely for them.")
}

# Predictor scaling check
# Check for NAs
if (sum(is.na(imsom_input_list$occ.covs)) > 0) stop("NAs found in occ.covs!")

# Check range (should be roughly -3 to +3)
summary(imsom_input_list$occ.covs)

# Replicate structure check
# Check detection covariates for Source 2 (Opportunistic)
str(imsom_input_list$det.covs[[2]])

# 6. Save
write_rds(imsom_input_list, here(proc_dir, "imsom_input_list.rds"))

if (!is.null(source1)) message(paste("Source 1 (ArHa) Sites:", dim(source1$y)[2], "| Replicates:", dim(source1$y)[3]))
if (!is.null(source2)) message(paste("Source 2 (Opp) Sites:", dim(source2$y)[2], "| Replicates:", dim(source2$y)[3]))
