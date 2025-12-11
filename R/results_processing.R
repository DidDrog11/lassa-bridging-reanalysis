################################################################################
## SCRIPT: results_processing.R
## PURPOSE: Generate spatial predictions (Dm), uncertainty maps, and 
##          coefficient plots from the fitted IMSOM.
################################################################################

# 1. Setup --------------------------------------------------------------------
source(here::here("packages.R"))

# Define directories
proc_dir <- here("data", "processed")
maps_dir <- here("results", "maps")
figs_dir <- here("results", "figures")
raw_dir  <- here("data", "raw")
if (!dir.exists(maps_dir)) dir.create(maps_dir, recursive = TRUE)
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# Load Model
model_path <- here("results", "models", "imsom_model_fit.rds")
if (!file.exists(model_path)) stop("Model file not found.")
out_imsom <- readRDS(model_path)

# Load Predictor Stack
stack_path <- here("data", "processed", "final_predictor_stack_scaled.tif")
pred_stack_scaled <- terra::rast(stack_path)

# Read in basinski raster to use as a mask
mask_path <- here("data", "raw", "basinski_comparison", "host_raster.tif")
basinski_aligned <- rast(mask_path)|> 
  terra::project(pred_stack_scaled, method = "near")

# Apply Mask
# This sets all pixels outside WA to NA in the predictor stack
pred_stack_masked <- terra::mask(pred_stack_scaled, basinski_aligned)

# 2. Prepare Prediction Data (Replicate Transformation) -----------------------
# Convert to Dataframe
# cells = TRUE: We keep the cell ID so we can map results back to the raster later
# na.rm = TRUE: This removes the ~100k NA pixels
pred_df <- terra::as.data.frame(pred_stack_masked, xy = TRUE, cells = TRUE, na.rm = TRUE)
# Occurrence formula
occ_formula <- ~ Tmu + Pmu + Nmu + Elev + 
  LC_12_Density +  # Croplands (Food resource)
  Pop +            # Linear Human Density (Commensalism)
  I(Pop^2)         # Quadratic Human Density (Exclusion at high density)
# Create Design Matrix (X.0)
X_0 <- model.matrix(occ_formula, data = pred_df)
coords_0 <- pred_df[, c("x", "y")]
valid_cell_ids <- pred_df$cell

# Load model input
imsom_input <- readRDS(here(proc_dir, "imsom_input_list.rds"))
species_names <- imsom_input$species_names

# 3. Run Spatial Prediction (Chunked) -----------------------------------------
output_data_path <- here(proc_dir, "prediction_results_multispecies.rds")

if(!file.exists(output_data_path)) {
  # A. Setup Chunks
  n_pixels <- nrow(X_0)
  n_species <- out_imsom$N
  chunk_size <- 5000
  n_chunks <- ceiling(n_pixels / chunk_size)
  
  # B. Initialize Output Vectors
  # Matrix: [Pixels x Species]
  psi_means_mat <- matrix(NA, nrow = n_pixels, ncol = n_species)
  colnames(psi_means_mat) <- species_names
  
  # Vector: [Pixels]
  target_species <- "Mastomys natalensis"
  sp_idx <- which(species_names == target_species)
  psi_sd_target <- numeric(n_pixels)
  
  # C. The Loop
  message("Starting Multi-Species Prediction Loop...")
  pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
  
  for (i in 1:n_chunks) {
    # Define indices
    start_idx <- (i - 1) * chunk_size + 1
    end_idx   <- min(i * chunk_size, n_pixels)
    current_idx <- start_idx:end_idx
    
    # Extract partial Design Matrix
    X_subset <- X_0[current_idx, , drop = FALSE]
    
    # Run Predict
    out_pred_chunk <- predict(out_imsom, X.0 = X_subset)
    
    # Calculate means for all species
    # shape: [samples, species, pixels]
    # Transpose to [pixels, species]
    chunk_means <- apply(out_pred_chunk$psi.0.samples, c(2, 3), mean)
    psi_means_mat[current_idx, ] <- t(chunk_means)
    
    # Calculate SD for Target Species Only
    # Extract target slice: [samples, 1, pixels]
    target_post <- out_pred_chunk$psi.0.samples[, sp_idx, ]
    if(length(current_idx) == 1) {
      psi_sd_target[current_idx] <- sd(target_post)
    } else {
      psi_sd_target[current_idx] <- apply(target_post, 2, sd)
    }
    
    rm(out_pred_chunk, target_post, chunk_means)
    gc()
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Save everything in one object
  saveRDS(list(valid_cells = valid_cell_ids,
               all_means = psi_means_mat,
               target_sd = psi_sd_target), output_data_path)
  
} else {
  message("Loading cached prediction results...")
}

# Load Data Back for Processing
pred_results <- readRDS(output_data_path)

# 4. Map back to Raster -------------------------------------------------------
r_template <- pred_stack_scaled[[1]]
terra::values(r_template) <- NA

# Loop through all species and save Mean Occupancy (Dm)
sp_names <- colnames(pred_results$all_means)

for(sp in sp_names) {
  clean_name <- gsub(" ", "_", sp)
  # Create Raster
  r_sp <- r_template
  r_sp[pred_results$valid_cells] <- pred_results$all_means[, sp]
  names(r_sp) <- paste0("Dm_", clean_name)
  # Save
  terra::writeRaster(r_sp, here(maps_dir, paste0("Dm_", clean_name, ".tif")), overwrite = TRUE)
}

# Uncertainty for Mastomys natalensis
r_unc <- r_template
r_unc[pred_results$valid_cells] <- pred_results$target_sd
names(r_unc) <- "Dm_Mastomys_natalensis_SD"
terra::writeRaster(r_unc, here(maps_dir, "Dm_Mastomys_natalensis_SD.tif"), overwrite = TRUE)

# 5. Visualisation: Coefficient Forest Plot -----------------------------------
message("Generating Coefficient Plot from Posterior Samples...")

# 1. Extract Beta Samples (Posterior Distribution)
# out_imsom$beta.samples is a matrix: [Samples x (Covariates * Species)]
beta_mat <- as.matrix(out_imsom$beta.samples)

# 2. Calculate Summary Statistics
beta_summary <- data.frame(Parameter = colnames(beta_mat),
                           Mean = apply(beta_mat, 2, mean),
                           Lower = apply(beta_mat, 2, quantile, probs = 0.025),
                           Upper = apply(beta_mat, 2, quantile, probs = 0.975)) |>
  tibble::rownames_to_column("row_id") |>
  # Filter out Intercepts to focus on environmental drivers
  filter(!grepl("\\(Intercept\\)", Parameter)) |>
  # Parse names: spOccupancy uses "Covariate-Species" format
  # We use extra="merge" to handle species names containing hyphens/spaces
  tidyr::separate(Parameter, into = c("Covariate", "Species"), sep = "-", extra = "merge")

# 3. Formatting for Plot
plot_data <- beta_summary |>
  mutate(Species = gsub("_", " ", Species), # Clean names (e.g., Mastomys_natalensis)
         # Define highlighting groups for the plot
         Highlight = case_when(Species == "Mastomys natalensis" ~ "Target",
                               Species %in% c("Rattus rattus", "Mus musculus") ~ "Invasive",
                               TRUE ~ "Native"))

# 4. Generate Plot
p_coef <- ggplot(plot_data, aes(x = Mean, y = reorder(Species, Mean), color = Highlight)) +
  # Zero line for significance
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # Whiskers (95% CI)
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 0.8) +
  # Dot (Mean)
  geom_point(size = 2.5) +
  # Facet by Covariate so we compare Species responses to the same driver
  facet_wrap(~Covariate, scales = "free_x", nrow = 2) +
  # Styling
  scale_color_manual(values = c("Target" = "#D55E00", "Invasive" = "#0072B2", "Native" = "grey60")) +
  theme_bw(base_size = 14) +
  labs(title = "Drivers of Rodent Occupancy in West Africa",
       subtitle = "Posterior Means and 95% Credible Intervals",
       x = "Effect Size (Logit Scale)",
       y = NULL,
       color = "Group") +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# Save
ggsave(here(figs_dir, "occupancy_coefficients.png"), p_coef, width = 14, height = 10)
message("Coefficient plot saved.")

# 6. Visualisation - Species Correlation Matrix (Spatial) -----------------
cor_mat_saved <- here(proc_dir, "cor_matrix.rds")

if (!file.exists(cor_mat_saved)) {
  message("Calculating correlation matrix from prediction results...")
  
  # Extract the matrix we calculated in the main loop
  mat_for_cor <- pred_results$all_means_mat
  
  # Calculate correlation
  cor_mat <- cor(mat_for_cor)
  saveRDS(cor_mat, cor_mat_saved)
} else {
  cor_mat <- readRDS(cor_mat_saved)
}

# Formatting
colnames(cor_mat) <- gsub("_", " ", colnames(cor_mat))
rownames(cor_mat) <- gsub("_", " ", rownames(cor_mat))

write.csv(cor_mat, here(figs_dir, "spatial_correlation_matrix.csv"))

# Plot Heatmap
if (require(corrplot)) {
  png(here(figs_dir, "species_correlation_heatmap.png"), width = 800, height = 800)
  corrplot::corrplot(cor_mat, 
                     method = "color", 
                     type = "lower", 
                     addCoef.col = "black", 
                     tl.col = "black", 
                     diag = FALSE, 
                     cl.lim = c(-1, 1), 
                     col = COL2('RdBu', 10), 
                     title = "Spatial Co-occurrence (Predicted)", 
                     mar = c(0, 0, 2, 0))
  dev.off()
}

