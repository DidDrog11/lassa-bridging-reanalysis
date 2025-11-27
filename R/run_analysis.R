################################################################################
## SCRIPT: run_analysis.R
##
## PURPOSE: Fit the Integrated Multi-Species Occupancy Model (IMSOM).
##          1. Load formatted inputs (imsom_input_list.rds).
##          2. Configure MCMC parameters (chains, iterations).
##          3. Run intMsPGOcc() from spOccupancy.
##          4. Save model object.
################################################################################

# 1. Setup --------------------------------------------------------------------
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

# Check Species Names (for reference)
print(imsom_input$species_names)