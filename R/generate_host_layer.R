################################################################################
## SCRIPT:  generate_host_layer.R
##
## PURPOSE: Fit the Integrated Multi-Species Occupancy Model (IMSOM).
##          1. Load formatted inputs (imsom_input_list.rds).
##          2. Configure MCMC parameters (chains, iterations).
##          3. Run intMsPGOcc() from spOccupancy.
##          4. Save model object.
################################################################################

# 1. Setup --------------------------------------------------------------------
# Load core packages (defined in packages.R)
source(here::here("packages.R"))

# Define directories
proc_dir <- here("data", "processed")
results_dir <- here("results", "models")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Load Inputs
imsom_input <- readRDS(here(proc_dir, "imsom_input_list.rds"))

# Extract components for clearer code
data_list <- list(
  y = imsom_input$y,
  occ.covs = imsom_input$occ.covs,
  det.covs = imsom_input$det.covs,
  sites = imsom_input$sites,
  coords = imsom_input$coords
)

# Species Names
spp_vector <- imsom_input$species_names # The vector of 8 names
data_list$species <- list(spp_vector, spp_vector) # List of length 2

# 2. Define Model Formulas and Priors ----------------------------------------------
# A. Ecological Process (Occupancy: psi)
# Strategy: 
# 1. Original Paper's Top Factors (Seasonality & Climate): Nm, Ncv, Pmu, Pcv, Tmu
# 2. Original Paper's Key Habitat: LC_9 (Savanna)
# 3. JSDM Hypothesis Drivers (Competition): Pop, LC_13 (Urban), LC_2 (Forest, although not a stable class so removed from predictors)

occ_formula <- ~ Tmu + Pmu + Pcv + Nm + Ncv + 
  LC_9_Density +  # Savanna (M. natalensis habitat)
  LC_13_Density + # Urban (Competition battleground)
  Pop             # Human Pressure (Commensal gradient)

# 1. Dropped 'LC_13_Density': Redundant with Pop and less informative for village settings.
# 2. Added 'I(Pop^2)': Quadratic term to capture the "hump-shaped" response 
#    (Low in Forest -> High in Village -> Low in City).

occ_formula_reduced <- ~ Tmu + Pmu + Nmu + Elev + 
  LC_12_Density +  # Croplands (Food resource)
  Pop +            # Linear Human Density (Commensalism)
  I(Pop^2)         # Quadratic Human Density (Exclusion at high density)

# B. Observation Process (Detection: p)
# We must define a separate formula for each data source in the list.

# Source 1: ArHa (High Quality)
# - effort: Continuous z-score of trap nights
# - quality: Binary flag (1 = Known, 0 = Imputed) to adjust intercept for uncertainty
det_formula_arha <- ~ effort + quality

# Source 2: Opportunistic (GBIF + Original)
# No effort data available (effort set to 0). We fit an intercept-only model.
det_formula_opp  <- ~ 1 

# Combine into a list
det_formulas <- list(det_formula_arha, det_formula_opp)

# C. Define Informative Priors (Detection)
# The ArHa dataset contains "Reporting Bias" (File Drawer Problem). 
# Targeted studies likely omitted low-effort/zero-catch nights.
# Therefore, Detection (p) is likely higher than in unbiased random surveys.
# We relax the priors to allow the model to find this high detection rate.

# B. Define Means (Neutral)
# Center at 0 (50% probability) instead of -3.42 (3% probability).
mu_alpha_arha <- c(0, 0, 0) # Intercept, Effort, Quality
mu_alpha_opp  <- c(-3.42)   # Keep Opportunistic low (it really is sparse)

# C. Define Variances (Wide/Weakly Informative)
# Increase variance to 4.0 (allows range roughly -4 to +4 on logit scale)
# This gives the data freedom to determine the slope and intercept.
var_alpha_arha <- c(4.0, 4.0, 4.0) 
var_alpha_opp  <- c(1)

# D. Occupancy Priors
# Previous: mean = 0, var = 2.72
# New: mean = -1.5, var = 1.0
# Rationale: Anchors the community near the observed average (-1.5), preventing
# Mastomys from drifting into the "Ubiquitous" (+2.0) space unchecked.
beta_comm_mean <- -1.5 
beta_comm_var  <- 1.0

# Structure per Documentation:
# alpha.comm.normal = list(mean = list(Source1, Source2), var = list(Source1, Source2))
priors_list <- list(
  alpha.comm.normal = list(
    mean = list(mu_alpha_arha, mu_alpha_opp),
    var  = list(var_alpha_arha, var_alpha_opp)
  ),
  beta.comm.normal = list(mean = beta_comm_mean, var = beta_comm_var),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig = list(a = list(0.1, 0.1), b = list(0.1, 0.1))
)

# 3. Configure MCMC & Initial Values ------------------------------------------
# A. Latent Factors (q)
# This determines the dimensionality of the residual correlation matrix.
# q=2 allows species to cluster along two "hidden" ecological axes (e.g., Forest vs. Commensal).
n_factors <- 2 

# B. MCMC Settings
n_samples    <- 50000
n_burn       <- 20000
n_thin       <- 10
n_chains     <- 3      

# C. Z Initialisation Values
get_max_obs <- function(y_arr) { apply(y_arr, c(1, 2), max, na.rm = TRUE) }

# Initialise z matrix (Species x Master Sites)
J_total <- nrow(data_list$occ.covs)
N_spp   <- length(spp_vector)
z_init  <- matrix(0, nrow = N_spp, ncol = J_total)

# Fill z based on Source 1 (ArHa)
sites_arha <- data_list$sites[[1]]
obs_arha   <- get_max_obs(data_list$y[[1]])
z_init[, sites_arha] <- pmax(z_init[, sites_arha], obs_arha, na.rm = TRUE)

# Fill z based on Source 2 (Opportunistic)
sites_opp <- data_list$sites[[2]]
obs_opp   <- get_max_obs(data_list$y[[2]])
z_init[, sites_opp] <- pmax(z_init[, sites_opp], obs_opp, na.rm = TRUE)

# Define Init List
inits_list <- list(
  alpha.comm = list(0, 0),   # Community mean detection (list of length n.data)
  beta.comm = 0,             # Community mean occupancy
  tau.sq.beta = 1,           # Community variance occupancy
  tau.sq.alpha = list(1, 1), # Community variance detection (list)
  beta = 0,                  # Species-specific occupancy params (start at mean)
  # Initialise alpha matrices: (Species x Params)
  # Source 1 has 3 params (Int, Effort, Qual), Source 2 has 1 (Int)
  alpha = list(matrix(0, nrow = N_spp, ncol = 3), matrix(0, nrow = N_spp, ncol = 1)),
  z = z_init)

# 4. Run Integrated JSDM ------------------------------------------------------
# Run the model
# NOTE: Using 'intMsPGOcc' (Latent Factor JSDM). 
# This model is structured based on the original manuscript, unfortunately there are convergence problems for the primary
# species of interest. We therefore use a reduced model which is more parsimonious.
# out_imsom <- intMsPGOcc(
#   occ.formula = occ_formula,
#   det.formula = det_formulas,
#   data = data_list,
#   inits = inits_list,
#   priors = priors_list,
#   n.samples = n_samples,
#   n.burn = n_burn,
#   n.thin = n_thin,
#   n.chains = n_chains,
#   n.omp.threads = 6,
#   verbose = TRUE,
#   n.report = 1000
# )

message("Starting model run (Host Layer)...")

out_imsom_reduced <- intMsPGOcc(
  occ.formula = occ_formula_reduced,
  det.formula = det_formulas,
  data = data_list,
  inits = inits_list,
  priors = priors_list,
  n.samples = n_samples,
  n.burn = n_burn,
  n.thin = n_thin,
  n.chains = n_chains,
  n.omp.threads = 6,
  verbose = TRUE,
  n.report = 1000
)

message("Model run complete")
# waicOcc(out_imsom)
# elpd          pD        WAIC 
# -1272.79986   99.16567  2743.93106 
waicOcc(out_imsom_reduced)
# elpd          pD        WAIC 
# -1231.34190   76.15738  2614.99857 

# 5. Save and Quick Diagnostics -----------------------------------------------
# A. Save Model Object
save_path <- here(results_dir, "imsom_model_fit.rds")
saveRDS(out_imsom_reduced, save_path)

# B. Convergence Check (R-hat)
# R-hat should be < 1.1 for convergence
cat("Beta (Occupancy):", max(out_imsom_reduced$rhat$beta, na.rm = TRUE), "\n")
cat("Alpha (Detection):", max(unlist(out_imsom_reduced$rhat$alpha), na.rm = TRUE), "\n")
