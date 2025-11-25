# --- QC CHECK: Final Predictor Stack ---

# 1. Calculate Global Statistics
# We verify Min, Max, and Mean to catch scaling/outlier errors
stats <- terra::global(Final_Predictor_Stack, fun = c("min", "max", "mean", "isNA"), na.rm = TRUE)
stats$Variable <- rownames(stats)
stats <- stats[, c("Variable", "min", "max", "mean", "isNA")] # Reorder

print(stats)

# --- 2. Define Passing Criteria (Mental Check) ---

# A. Climate Means (Tmu, Pmu)
# Tmu: Should be approx 15 to 35 (Degrees C).
# Pmu: Should be approx 0 to 20 (mm/day). 

# B. Indices (Nmu, Ncv, Pcv)
# Nmu: Must be -0.2 to 1.0.
# CVs (Pcv, Ncv): Must be > 0. High values (>1.0) are possible in deserts but shouldn't be negative.

# C. Colwell's & Probabilities (Pc, Pm, Densities)
# Pc, Pm, LC_Densities: MUST be strictly between 0 and 1.
# If Max > 1.00001 -> Failure.

# D. Duration
# Pdur, Ndur: Must be between 0 and 12 (months).

# --- 3. Visual Check of Gradients ---

# Set up plotting area
par(mfrow = c(2, 2))

# Check 1: The North-South Gradient (NDVI vs Precip)
terra::plot(Final_Predictor_Stack[["Nmu"]], main = "Mean NDVI (Check N-S Gradient)")
terra::plot(Final_Predictor_Stack[["Pmu"]], main = "Mean Precip (Check N-S Gradient)")

# Check 2: Seasonality Logic
# Contingency (Pm) should be high in the Savanna belt, low in deep desert/forest
terra::plot(Final_Predictor_Stack[["Pm"]], main = "Precip Seasonality (Pm)")

# Check 3: Urban Density (Check for localized hotspots)
if ("LC_13_Density" %in% names(Final_Predictor_Stack)) {
  terra::plot(Final_Predictor_Stack[["LC_13_Density"]], main = "Urban Density")
} else {
  plot(Final_Predictor_Stack[["Pop"]], main = "Pop Density")
}

par(mfrow = c(1, 1))

# --- 4. Correlation Check ---
# Ecological: NDVI should correlate positively with Precipitation in W. Africa
# Extract a sample of 5000 points to check correlation quickly
samp <- terra::spatSample(Final_Predictor_Stack, size = 5000, na.rm = TRUE, as.df = TRUE)

cor_np <- cor(samp$Nmu, samp$Pmu)
message(paste("Correlation between Nmu and Pmu:", round(cor_np, 3)))
