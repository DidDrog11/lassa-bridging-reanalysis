################################################################################
## SCRIPT: generate_pathogen_layer.R
## PURPOSE: Fit Boosted Regression Tree (BRT) for Pathogen Risk (Dl).
##          Model: Pr(LASV+ | Mastomys Presence) ~ Environment
## INPUTS:  data/processed/pathogen_training_data.rds (Cleaned & Blocked)
##          data/processed/final_predictor_stack_scaled.tif
## OUTPUTS: results/maps/Dl_Pathogen_Layer.tif
################################################################################

# 1. Setup --------------------------------------------------------------------
source(here::here("packages.R"))

# Define directories
proc_dir <- here("data", "processed")
maps_dir <- here("results", "maps")
models_dir <- here("results", "models")
raw_dir <- here("data", "raw")
if (!dir.exists(maps_dir)) dir.create(maps_dir, recursive = TRUE)
if (!dir.exists(models_dir)) dir.create(models_dir, recursive = TRUE)

# Load Training Data 
# This was created in clean_data.R Section C
# It already has: outcome (0/1), spatial blocking applied, and predictors extracted
train_data <- readRDS(here(proc_dir, "pathogen_training_data.rds"))

# Load the scaled Predictor Stack
stack_path <- here(proc_dir, "final_predictor_stack_scaled.tif")
if (!file.exists(stack_path)) stop("Scaled predictor stack not found.")
pred_stack <- terra::rast(stack_path)

# Assemble ecological predictor stack
# Produced from the IMSOM extracted in results_processing.R
r_rattus_path <- here(maps_dir, "Dm_Rattus_rattus.tif")
r_mus_path    <- here(maps_dir, "Dm_Mus_musculus.tif")
if (!file.exists(r_rattus_path) || !file.exists(r_mus_path)) {
  stop("Invasive species maps not found. Run results_processing.R first.")
}
r_rattus <- terra::rast(r_rattus_path)
r_mus    <- terra::rast(r_mus_path)
r_rattus <- terra::resample(r_rattus, pred_stack, method = "near")
r_mus    <- terra::resample(r_mus, pred_stack, method = "near")

# Create final DL stack
Dl_stack <- c(pred_stack, r_rattus, r_mus)
# Ensure names are clean for the BRT formulas
names(Dl_stack) <- c(names(pred_stack), "Dm_Rattus", "Dm_Mus")

# Re-extract the training data now that the ecological predictors are added
train_points_raw <- train_data |>
  dplyr::select(latitude, longitude, n_pos, n_tested, outcome, source) |>
  filter(!is.na(latitude))

v_pts <- terra::vect(train_points_raw, geom = c("longitude", "latitude"), crs = "EPSG:4326")
extract_vals <- terra::extract(Dl_stack, v_pts, ID = FALSE)

train_data <- cbind(train_points_raw, extract_vals) |>
  na.omit()

# 2. Configure BRT Model ------------------------------------------------------
# Define Predictors
metadata_cols <- c("latitude", "longitude", "n_tested", "n_pos", "outcome", "source")
pred_vars <- c(
  # 1. Climate/Seasonality (Basinski Baseline)
  "Pm", "Pc", "Pmax", "Pmu", "Tmu",
  
  # 2. Vegetation
  "Nmax", "Nmin", "Nmu", "Ndur",
  
  # 3. Habitat (Land Cover)
  "LC_10_Density", "LC_12_Density", "LC_9_Density",
  
  # 4. The Socio-economic Shield (Infrastructure)
  "NTL", 
  
  # 5. The Biotic Shield (Competitors)
  "Dm_Rattus", "Dm_Mus"
)
# Check if all vars exist in data
missing_vars <- setdiff(pred_vars, names(train_data))
if(length(missing_vars) > 0) stop(paste("Missing variables:", paste(missing_vars, collapse=", ")))

# 3. Fit Boosted Regression Tree (BRT) ----------------------------------------
# Set seed for reproducibility
set.seed(123)

# A. Fit Single "Best" Model (Full Data)
# This is used for Variable Importance and the Primary Map
# Learning Rate Logic: If N < 50, use slow learning rate
lr_main <- if (nrow(train_data) < 50) 0.0005 else 0.001

brt_full <- dismo::gbm.step(
  data = as.data.frame(train_data),
  gbm.x = which(names(train_data) %in% pred_vars),
  gbm.y = which(names(train_data) == "outcome"),
  family = "bernoulli",
  tree.complexity = 3,
  learning.rate = lr_main,
  bag.fraction = 0.75,
  n.folds = 10,
  verbose = TRUE
)

if (!is.null(brt_full)) {
  saveRDS(brt_full, here(models_dir, "pathogen_model_fit_full.rds"))
}

# B. Fit Bootstrap Models (Uncertainty Estimation)
n_boots <- 25
lr_boot <- 0.001 

# Initialize accumulators using the correct stack
sum_pred <- terra::rast(Dl_stack[[1]]); terra::values(sum_pred) <- 0
sum_sq_pred <- terra::rast(Dl_stack[[1]]); terra::values(sum_sq_pred) <- 0

message("Starting Bootstraps...")

for (i in 1:n_boots) {
  message(paste("Bootstrap iteration:", i, "/", n_boots))
  
  # Resample
  boot_idx <- sample(nrow(train_data), nrow(train_data), replace = TRUE)
  boot_data <- train_data[boot_idx, ]
  
  # Fit Model
  # Use tryCatch to skip folds that fail to converge (common with small N)
  brt_boot <- tryCatch({
    dismo::gbm.step(
      data = as.data.frame(boot_data),
      gbm.x = which(names(boot_data) %in% pred_vars),
      gbm.y = which(names(boot_data) == "outcome"),
      family = "bernoulli",
      tree.complexity = 2,
      learning.rate = lr_boot, 
      bag.fraction = 0.75,
      n.folds = 10,
      verbose = FALSE, 
      plot.main = FALSE
    )
  }, error = function(e) return(NULL)) 
  
  if (!is.null(brt_boot)) {
    # Predict using Dl_stack (NOT pred_stack)
    pred_r <- terra::predict(Dl_stack, brt_boot, 
                             n.trees = brt_boot$gbm.call$best.trees, 
                             type = "response", na.rm=TRUE)
    
    # Accumulate
    sum_pred <- sum_pred + pred_r
    sum_sq_pred <- sum_sq_pred + (pred_r^2)
  } else {
    message("  - Bootstrap failed/skipped.")
  }
}

# 4. Generate Final Rasters (Mean & Uncertainty) ------------------------------
# A. Primary Map (Best Model)
Dl_primary <- terra::predict(Dl_stack, brt_full, 
                             n.trees = brt_full$gbm.call$best.trees, 
                             type = "response", na.rm = TRUE)
names(Dl_primary) <- "Dl_Pathogen_Risk_Primary"

# B. Ensemble Mean
Dl_ensemble_mean <- sum_pred / n_boots
names(Dl_ensemble_mean) <- "Dl_Pathogen_Risk_Ensemble"

# C. Uncertainty (SD)
var_rast <- (sum_sq_pred - (sum_pred^2 / n_boots)) / (n_boots - 1)
Dl_sd <- sqrt(var_rast)
names(Dl_sd) <- "Dl_Pathogen_Uncertainty"

# Save
terra::writeRaster(Dl_primary, here(maps_dir, "Dl_Pathogen_Layer_Primary.tif"), overwrite = TRUE)
terra::writeRaster(Dl_ensemble_mean, here(maps_dir, "Dl_Pathogen_Layer_Ensemble.tif"), overwrite = TRUE)
terra::writeRaster(Dl_sd, here(maps_dir, "Dl_Pathogen_Uncertainty.tif"), overwrite = TRUE)