################################################################################
## SCRIPT: calculate_risk_and_validate_comparison.R
## PURPOSE: 1. Calculate Composite Risk (Dx) for scenarios (Custom vs Basinski).
##          2. Apply IUCN Range Mask (Sanity Filter).
##          3. Calibrate against Human Seroprevalence (GLM).
##          4. Estimate Annual Incidence (SIRS) with Lambda Sensitivity.
##          5. Validate and Compare Performance.
################################################################################

# 1. Setup --------------------------------------------------------------------
source(here::here("packages.R"))
# Directories
maps_dir <- here("results", "maps")
data_dir <- here("data", "raw")
proc_dir <- here("data", "processed")
figs_dir <- here("results", "figures")
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

Pop_count <- terra::rast(here(proc_dir, "pop_count_2020_05deg.tif"))


# 2. Inputs for Comparison ------------------------------------------------
# A. Load IUCN Range Mask
# This acts as the final "Biological Reality Check"
iucn_path <- here(data_dir, "iucn", "data_0.shp")
if(file.exists(iucn_path)) {
  iucn_vect <- terra::vect(iucn_path)
  # Rasterize to match Population Grid (1 = Inside Range, NA = Outside)
  iucn_mask <- terra::rasterize(iucn_vect, Pop_count, field = 1, background=NA)
} else {
  warning("IUCN Shapefile not found. Proceeding without range mask.")
  iucn_mask <- Pop_count # Dummy
  values(iucn_mask) <- 1
}

# B. Custom Layers
dm_custom_path <- here(maps_dir, "Dm_Mastomys_natalensis.tif")
dl_custom_path <- here(maps_dir, "Dl_Pathogen_Layer_Ensemble.tif")
if(!file.exists(dl_custom_path)) dl_custom_path <- here(maps_dir, "Dl_Pathogen_Layer.tif")

# C. Basinski Layers
dm_basinski_path <- here(data_dir, "basinski_comparison", "host_raster.tif")
dl_basinski_path <- here(data_dir, "basinski_comparison", "pathogen_raster.tif")

# D. Generate Death Rate Raster (drast)
# Approx Life Expectancy 2023 (World Bank approx) -> d = 1/LE
le_data <- data.frame(GID_0 = c("BEN", "BFA", "CIV", "GHA", "GIN", "GMB",
                                "GNB", "LBR", "MLI", "MRT", "NER", 
                                "NGA", "SEN", "SLE", "TGO"),
                      LE = c(61, 61, 62, 65, 61, 63, 64, 62, 60, 68, 61, 54, 69, 62, 63))

# Fetch boundaries and join LE
wa_adm0_vect <- gadm(country = le_data$GID_0, level = 0, path = data_dir)
# Match IDs safely
values(wa_adm0_vect) <- values(wa_adm0_vect) %>% 
  left_join(le_data, by="GID_0") %>%
  mutate(DeathRate = 1/LE)

drast <- terra::rasterize(wa_adm0_vect, Pop_count, field = "DeathRate", background = NA)

# E. Define Scenarios
scenarios <- list()

if(file.exists(dm_custom_path) & file.exists(dl_custom_path)) {
  scenarios[["Full_Custom"]] <- list(dm = dm_custom_path, dl = dl_custom_path)
}
if(file.exists(dm_custom_path) & file.exists(dl_basinski_path)) {
  scenarios[["Hybrid_OurDm_BasDl"]] <- list(dm = dm_custom_path, dl = dl_basinski_path)
}
if(file.exists(dm_basinski_path) & file.exists(dl_custom_path)) {
  scenarios[["Hybrid_BasDm_OurDl"]] <- list(dm = dm_basinski_path, dl = dl_custom_path)
}
if(file.exists(dm_basinski_path) & file.exists(dl_basinski_path)) {
  scenarios[["Full_Basinski"]] <- list(dm = dm_basinski_path, dl = dl_basinski_path)
}

# F. Validation Data
human_data <- read.csv(here(proc_dir, "human_calibration_data.csv")) |>
  filter(!is.na(lat) & !is.na(lon) & n_test > 0)
human_vect <- terra::vect(human_data, geom = c("lon", "lat"), crs = "EPSG:4326")

moore_data <- tryCatch(read.csv(here(proc_dir, "validation_cases_admin2.csv")), error = function(e) NULL)
lga_status <- tryCatch(read.csv(here(proc_dir, "validation_binary_lga.csv")), error = function(e) NULL)

# 3. Comparison Loop ----------------------------------------------------
results_summary <- data.frame()
model_fits <- list()

for (scen_name in names(scenarios)) {
  message(paste0("\nProcessing: ", scen_name))
  
  # A. Load Layers
  paths <- scenarios[[scen_name]]
  Dm_raw <- terra::rast(paths$dm)
  Dl_raw <- terra::rast(paths$dl)
  
  # B. Alignment (Resample to Pop Grid)
  if (!terra::compareGeom(Pop_count, Dm_raw, stopOnError = FALSE)) {
    Dm_raster <- terra::resample(Dm_raw, Pop_count, method = "bilinear")
  } else { Dm_raster <- Dm_raw }
  
  if (!terra::compareGeom(Pop_count, Dl_raw, stopOnError = FALSE)) {
    Dl_raster <- terra::resample(Dl_raw, Pop_count, method = "bilinear")
  } else { Dl_raster <- Dl_raw }
  
  # C. Masking (Apply IUCN Range)
  # This sets risk to NA outside the known host range
  Dm_raster <- terra::mask(Dm_raster, iucn_mask)
  Dl_raster <- terra::mask(Dl_raster, iucn_mask)
  
  # D. Composite Risk (Dx)
  Dx_raster <- Dm_raster * Dl_raster
  names(Dx_raster) <- "Dx" 
  terra::writeRaster(Dx_raster, here(maps_dir, paste0("Dx_", scen_name, ".tif")), overwrite = TRUE)
  
  # E. Calibration: Testing the Socio-economic Shield (NTL)
  # -------------------------------------------------------
  calib_data <- human_data
  # Extract Risk (Dx)
  calib_data$Dx <- terra::extract(Dx_raster, human_vect, ID = FALSE)[,1]
  
  # Extract Shield (NTL)
  if(!exists("NTL_raster")) NTL_raster <- terra::rast(here(proc_dir, "final_predictor_stack_scaled.tif"))[["NTL"]]
  calib_data$NTL <- terra::extract(NTL_raster, human_vect, ID = FALSE)[,1]
  
  calib_data <- na.omit(calib_data |>
                          select(-id))
  calib_data$PropAb <- calib_data$n_pos / calib_data$n_test
  
  # Model 1: Quasi-Binomial (Base)
  glm_base <- glm(PropAb ~ Dx, family = quasibinomial(link = "logit"), 
                  weights = n_test, data = calib_data)
  
  # Model 2: Quasi-Binomial (Shielded)
  glm_shield <- glm(PropAb ~ Dx + NTL, family = quasibinomial(link = "logit"), 
                    weights = n_test, data = calib_data)
  
  # Prediction for Base (Unshielded)
  # (Risk ~ Climate only)
  Omega_Base <- terra::predict(Dx_raster, glm_base, type = "response", na.rm = TRUE)
  
  # Prediction for Shielded
  # (Risk ~ Ecology + Infrastructure)
  pred_stack_shield <- c(Dx_raster, NTL_raster)
  names(pred_stack_shield) <- c("Dx", "NTL")
  Omega_Shield <- terra::predict(pred_stack_shield, glm_shield, type = "response", na.rm = TRUE)
  
  # Statistical Test
  model_anova <- anova(glm_base, glm_shield, test = "F")
  p_val_ntl <- model_anova$`Pr(>F)`[2]
  ntl_coef <- coef(glm_shield)["NTL"]
  
  message(paste("  -> NTL P-value:", signif(p_val_ntl, 4)))
  
  # Decision Logic
  if (!is.na(p_val_ntl) && p_val_ntl < 0.05 && ntl_coef < 0) {
    message("  -> Shield Model Selected (Significant Protective Effect).")
    best_glm <- glm_shield
    model_type <- "Shielded"
    
    # Predict using Stack (Dx + NTL)
    pred_stack_calib <- c(Dx_raster, NTL_raster)
    names(pred_stack_calib) <- c("Dx", "NTL")
    Omega_raster <- terra::predict(pred_stack_calib, best_glm, type = "response", na.rm = TRUE)
    
  } else {
    message("  -> Base Model Selected (NTL not significant).")
    best_glm <- glm_base
    model_type <- "Base"
    Omega_raster <- terra::predict(Dx_raster, best_glm, type = "response", na.rm = TRUE)
  }
  
  # Store stats for output
  model_fits[[scen_name]] <- best_glm
  dev_expl <- 1 - (best_glm$deviance / best_glm$null.deviance)
  
  # Post-Processing & Incidence
  # Mask Seroprevalence by IUCN
  Omega_raster <- terra::mask(Omega_raster, iucn_mask)
  names(Omega_raster) <- paste0("Omega_", scen_name)
  terra::writeRaster(Omega_raster, here(maps_dir, paste0("Omega_", scen_name, ".tif")), overwrite = TRUE)
  
  # Incidence Calculation (SIRS)
  # n.b. mu is mislabelled here, it actually refers to f the infection specific fatality rate
  gam_val <- 12; mu_val <- 0.02
  Base_Cases <- Omega_raster * Pop_count * (drast + gam_val) * (drast) / (gam_val * (1 - mu_val))
  
  lambda_vals <- c(0, 0.03, 0.064)
  
  for (lam in lambda_vals) {
    
    # Common Multiplier
    mult_factor <- (drast + gam_val) * (drast + lam) / (gam_val * (1 - mu_val))
    
    # A. Calculate Base Incidence
    Inc_Base <- Omega_Base * Pop_count * mult_factor
    Inc_Base <- terra::mask(Inc_Base, iucn_mask)
    total_base <- terra::global(Inc_Base, "sum", na.rm=TRUE)$sum
    
    # B. Calculate Shielded Incidence
    Inc_Shield <- Omega_Shield * Pop_count * mult_factor
    Inc_Shield <- terra::mask(Inc_Shield, iucn_mask)
    total_shield <- terra::global(Inc_Shield, "sum", na.rm=TRUE)$sum
    
    # Store Results
    results_summary <- rbind(results_summary, data.frame(
      Scenario = scen_name, 
      Lambda = lam, 
      Model_Type = "Unshielded",
      Total_Infections = round(total_base, 0)
    ))
    
    results_summary <- rbind(results_summary, data.frame(
      Scenario = scen_name, 
      Lambda = lam, 
      Model_Type = "Shielded",
      Total_Infections = round(total_shield, 0)
    ))
  }
}

print(results_summary)
saveRDS(results_summary, here("results", "tables", "method_comparison_summary.rds"))

# Define NTL Contexts to visualize the "Shield" effect
# We assume:
# - Rural: NTL approx -0.13 (from your raw data min)
# - Urban Core: NTL approx 20 (from your synthetic anchors)
ntl_contexts <- data.frame(
  Context = c("Rural (Dark)", "Urban Core (Bright)"),
  NTL = c(-0.13, 20)
)

pred_frame_base <- data.frame(Dx = seq(0, 1, length.out = 100))
plot_lines <- data.frame()

for(scen in names(model_fits)) {
  model <- model_fits[[scen]]
  
  # Check if this model is "Shielded" (has NTL coefficient)
  has_shield <- "NTL" %in% names(coef(model))
  
  if(has_shield) {
    # If Shielded, we predict TWO curves to show the interaction
    for(i in 1:nrow(ntl_contexts)) {
      # Create prediction frame with specific NTL
      tmp_pred <- pred_frame_base
      tmp_pred$NTL <- ntl_contexts$NTL[i]
      
      tmp_out <- data.frame(
        Dx = tmp_pred$Dx,
        Predicted = predict(model, newdata = tmp_pred, type = "response"),
        Scenario = scen,
        Type = ntl_contexts$Context[i]
      )
      plot_lines <- rbind(plot_lines, tmp_out)
    }
    
  } else {
    # If Base Model (e.g., Basinski or if Shield was rejected), 
    # it predicts the same line regardless of NTL.
    # We map it to "Base Assumption"
    tmp_out <- data.frame(
      Dx = pred_frame_base$Dx,
      Predicted = predict(model, newdata = pred_frame_base, type = "response"),
      Scenario = scen,
      Type = "Base Assumption (No Shield)"
    )
    plot_lines <- rbind(plot_lines, tmp_out)
  }
}

# Generate the Plot
# We use Linetype to distinguish Rural vs Urban, and Color for Scenario
p_comp <- ggplot(plot_lines, aes(x = Dx, y = Predicted, color = Scenario)) +
  geom_line(aes(linetype = Type), size = 1.2) +
  scale_color_brewer(palette = "Set1") +
  scale_linetype_manual(values = c(
    "Rural (Dark)" = "solid", 
    "Urban Core (Bright)" = "dotted", 
    "Base Assumption (No Shield)" = "dashed"
  )) +
  theme_bw() +
  labs(title = "The Socio-economic Shield Effect",
       subtitle = "Calibration curves showing how Infrastructure (NTL) dampens Risk (Dx)",
       x = "Composite Ecological Hazard (Dx)",
       y = "Predicted Human Seroprevalence",
       linetype = "Context",
       color = "Model Scenario") +
  theme(legend.position = "bottom", legend.box = "vertical")

# Save
ggsave(here(figs_dir, "comparison_calibration_curves_shielded.png"), p_comp, width = 9, height = 7)
