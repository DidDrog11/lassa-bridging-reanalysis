################################################################################
## SCRIPT: sensitivity_analyses.R
## PURPOSE: Propagate uncertainty from Host (IMSOM) and Pathogen (BRT) layers 
##          to the final Risk (Dx) layer and attribute sources of error.
## INPUTS:  Dm_Mastomys_Mean.tif / Dm_Mastomys_SD.tif
##          Dl_Pathogen_Layer_Ensemble.tif / Dl_Pathogen_Uncertainty.tif
## OUTPUTS: Dx_Uncertainty.tif, Uncertainty_Source_Ratio.tif
################################################################################

source(here::here("packages.R"))

# Directories
maps_dir <- here("results", "maps")
figs_dir <- here("results", "figures")
data_dir <- here("data", "raw")
proc_dir <- here("data", "processed")
tabs_dir <- here("results", "tables")

# 1. Load Mean and Uncertainty Rasters ----------------------------------
# A. Host Layer (Mastomys)
# Mean Prediction
Dm_mean <- terra::rast(here(maps_dir, "Dm_Mastomys_natalensis.tif"))
# Uncertainty (SD) - Created in results_processing.R
Dm_sd   <- terra::rast(here(maps_dir, "Dm_Mastomys_natalensis_SD.tif"))

# B. Pathogen Layer (Lassa)
# Mean Prediction
Dl_mean <- terra::rast(here(maps_dir, "Dl_Pathogen_Layer_Ensemble.tif"))
# Uncertainty (SD) - Created in generate_pathogen_layer.R
Dl_sd   <- terra::rast(here(maps_dir, "Dl_Pathogen_Uncertainty.tif"))

# C. Align Rasters (Safety Check)
# Ensure extents match exactly before math
if (!terra::compareGeom(Dm_mean, Dl_mean, stopOnError = FALSE)) {
  Dl_mean <- terra::resample(Dl_mean, Dm_mean, method = "bilinear")
  Dl_sd   <- terra::resample(Dl_sd, Dm_mean, method = "bilinear")
}

# 2. Propagate Error to Composite Risk (Dx) -----------------------------
# Calculate Global Correlation (Pearson's r)
# When Host prob is high, is Pathogen prob high?
global_cor <- terra::layerCor(c(Dm_mean, Dl_mean), fun = "pearson", na.rm = TRUE)
print(global_cor$correlation)

# Local Correlation Map (Moving Window)
# Calculate correlation in a moving window (e.g., ~50km radius)
r_cor_local <- terra::focalPairs(c(Dm_mean, Dl_mean), w = 9, fun = "pearson", na.rm = TRUE)
names(r_cor_local) <- "Local_Correlation_Dm_Dl"

p_cor <- ggplot() +
  geom_spatraster(data = r_cor_local) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limit = c(-1, 1),
                       name = "Correlation (r)") +
  theme_minimal() +
  labs(title = "Local Correlation between Host & Pathogen",
       subtitle = "Red areas (High +ve r) imply errors might compound.\nBlue/White areas support independence.")

# Formula for Variance of a Product (assuming independence):
# Var(X*Y) = Var(X)*E[Y]^2 + Var(Y)*E[X]^2 + Var(X)*Var(Y)
# Where X = Dm, Y = Dl, and Var = SD^2
Var_Dm <- Dm_sd^2
Var_Dl <- Dl_sd^2

# Calculate Variance of Dx
Var_Dx <- (Var_Dm * Dl_mean^2) + (Var_Dl * Dm_mean^2) + (Var_Dm * Var_Dl)

# Convert back to SD
Dx_sd <- sqrt(Var_Dx)
names(Dx_sd) <- "Dx_SD"

terra::writeRaster(Dx_sd, here(maps_dir, "Dx_Uncertainty.tif"), overwrite = TRUE)

# Calculate Coefficient of Variation (CV) for context
# CV = SD / Mean. 
Dx_mean <- Dm_mean * Dl_mean
Dx_cv <- Dx_sd / Dx_mean
names(Dx_cv) <- "Dx_CV"
terra::writeRaster(Dx_cv, here(maps_dir, "Dx_CV.tif"), overwrite = TRUE)

# 3. Where does the error come from? ------------------------------------
# We calculate the contribution of each term to the total variance
# Term 1: Host Contribution ~ Var(Dm) * Dl^2
# Term 2: Pathogen Contribution ~ Var(Dl) * Dm^2
Contrib_Host <- (Var_Dm * Dl_mean^2)
Contrib_Pathogen <- (Var_Dl * Dm_mean^2)

# Create a Ratio Map
# > 1: Host is the main source of error
# < 1: Pathogen is the main source of error
# log-transform for symmetric plotting: 0 = Equal, +ve = Host, -ve = Pathogen
Uncertainty_Ratio <- log(Contrib_Host / Contrib_Pathogen)
names(Uncertainty_Ratio) <- "Log_Uncertainty_Ratio_Host_vs_Pathogen"

terra::writeRaster(Uncertainty_Ratio, here(maps_dir, "Uncertainty_Source_Ratio.tif"), overwrite = TRUE)

# 4. Visualisation ------------------------------------------------------
# Plot 1: Total Absolute Uncertainty (Where are we unsure?)
p_unc <- ggplot() +
  geom_spatraster(data = Dx_sd) +
  scale_fill_viridis_c(option = "inferno", name = "Risk SD", na.value = "transparent") +
  theme_minimal() +
  labs(title = "Absolute Uncertainty in Composite Risk (SD)",
       subtitle = "High values indicate regions where risk estimates are unstable")

ggsave(here(figs_dir, "Map_Dx_Uncertainty.png"), p_unc, width = 8, height = 6)

# Plot 2: Source of Uncertainty (Why are we unsure?)
p_source <- ggplot() +
  geom_spatraster(data = Uncertainty_Ratio) +
  scale_fill_gradient2(
    low = "#0072B2",   # Blue = Pathogen Uncertainty Dominates
    mid = "white",     # White = Equal
    high = "#D55E00",  # Red = Host Uncertainty Dominates
    midpoint = 0,
    name = "Log Ratio",
    na.value = "transparent"
  ) +
  theme_minimal() +
  labs(title = "Source of Risk Uncertainty",
       subtitle = "Red: Unsure about Rats (Host) | Blue: Unsure about Virus (Pathogen)")

ggsave(here(figs_dir, "Map_Uncertainty_Source.png"), p_source, width = 8, height = 6)


# 5. Monte Carlo Simulation: Incidence Uncertainty ------------------------
# A. Load Ancillary Data
# Population
Pop_count <- terra::rast(here(proc_dir, "pop_count_2020_05deg.tif"))

# NTL
NTL_raster <- terra::rast(here(proc_dir, "final_predictor_stack_scaled.tif"))[["NTL"]]

# Death Rate Raster
le_data <- data.frame(GID_0 = c("BEN", "BFA", "CIV", "GHA", "GIN", "GMB",
                                "GNB", "LBR", "MLI", "MRT", "NER", 
                                "NGA", "SEN", "SLE", "TGO"),
                      LE = c(61, 61, 62, 65, 61, 63, 64, 62, 60, 68, 61, 54, 69, 62, 63))

# Fetch boundaries and join LE
wa_adm0_vect <- gadm(country = le_data$GID_0, level = 0, path = data_dir)
# Match IDs safely
values(wa_adm0_vect) <- values(wa_adm0_vect) |>
  left_join(le_data, by = "GID_0") |>
  mutate(DeathRate = 1/LE)

drast <- terra::rasterize(wa_adm0_vect, Pop_count, field = "DeathRate", background = NA)

# Admin Zones for Aggregation
adm0_rast <- terra::rasterize(wa_adm0_vect, Pop_count, field = "COUNTRY")
adm0_names <- data.frame(ID = 1:length(unique(wa_adm0_vect$COUNTRY)), Name = sort(unique(wa_adm0_vect$COUNTRY)))

# Refit the Shielded GLM
human_aug <- read.csv(here(proc_dir, "human_calibration_data.csv")) |>
  filter(!is.na(lat) & !is.na(lon) & n_test > 0) |>
  select(-id)
human_vect <- terra::vect(human_aug, geom = c("lon", "lat"), crs = "EPSG:4326")

calib_data <- human_aug
calib_data$Dx <- terra::extract(Dx_mean, human_vect, ID = FALSE)[,1]
calib_data$NTL <- terra::extract(NTL_raster, human_vect, ID = FALSE)[,1]
calib_data <- na.omit(calib_data)

# Shielded Model
glm_shield <- glm(cbind(n_pos, n_test - n_pos) ~ Dx + NTL, 
                  family = quasibinomial(link = "logit"), data = calib_data)

# C. Simulation Loop
n_sims <- 100 
lambda_select <- 0.03
results_adm0 <- list()

# Load IUCN Mask
iucn_path <- here("data", "raw", "iucn", "data_0.shp")
if(file.exists(iucn_path)) {
  iucn_vect <- terra::vect(iucn_path)
  iucn_mask <- terra::rasterize(iucn_vect, Pop_count, field=1, background=NA)
} else {
  stop("IUCN Shapefile not found!")
}

pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

for(i in 1:n_sims) {
  
  # 1. Simulate Risk Surface (Dx)
  # Sample pixel from Normal(Mean, SD), clamped to [0,1]
  Dx_sim <- terra::app(c(Dx_mean, Dx_sd), fun = function(x) {
    val <- rnorm(1, mean = x[1], sd = x[2])
    return(max(0, min(1, val)))
  })
  names(Dx_sim) <- "Dx"
  
  # 2. Predict Prevalence (Omega)
  pred_stack <- c(Dx_sim, NTL_raster)
  names(pred_stack) <- c("Dx", "NTL")
  Omega_sim <- terra::predict(pred_stack, glm_shield, type = "response", na.rm = TRUE)
  Omega_sim <- terra::mask(Omega_sim, iucn_mask)
  
  # 3. Calculate Incidence
  gam_val <- 12; mu_val <- 0.02
  
  # Step A: Base Calculation
  Base_Cases <- Omega_sim * Pop_count * (drast + gam_val) * (drast) / (gam_val * (1 - mu_val))
  
  # Step B: Lambda Multiplier (Includes / drast)
  # Note: The drast in Base and the drast denominator here cancel out algebraically,
  lam_mult_rast <- (drast + lambda_select) / drast
  
  Incidence_sim <- Base_Cases * lam_mult_rast
  
  # 4. Aggregate by Country
  z_adm0 <- terra::zonal(Incidence_sim, adm0_rast, fun = "sum", na.rm = TRUE)
  z_adm0$Sim_ID <- i
  results_adm0[[i]] <- z_adm0
  
  setTxtProgressBar(pb, i)
}
close(pb)

df_adm0_all <- do.call(rbind, results_adm0)
colnames(df_adm0_all)[1] <- "ID"

summary_adm0 <- df_adm0_all |>
  rename(Cases = 2) |>
  group_by(ID) |>
  summarise(Mean_Cases = round(mean(Cases, na.rm = TRUE)),
            Lower_CI = round(quantile(Cases, 0.025, na.rm = TRUE)),
            Upper_CI = round(quantile(Cases, 0.975, na.rm=TRUE)),
            SD_Cases = round(sd(Cases, na.rm=TRUE))) |>
  mutate(Country = unique(df_adm0_all$ID)) |>
  select(Country, Mean_Cases, Lower_CI, Upper_CI, SD_Cases) |>
  arrange(desc(Mean_Cases))

write.csv(summary_adm0, here(tabs_dir, "incidence_uncertainty_country.csv"), row.names = FALSE)

print(head(summary_adm0))