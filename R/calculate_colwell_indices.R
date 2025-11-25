# R/calculate_colwell_indices.R

# Function to calculate Colwell's Indices (Constancy and Contingency) 
# for a single pixel's time series.
# This function is designed to be called by terra::app().

#' Calculates Colwell's Indices (Pc and Pm)
#'
#' @param x A vector of chronological monthly values (e.g., precipitation or NDVI).
#' @return A numeric vector of length 2: c(Constancy (Pc), Contingency (Pm)).
#' @references Colwell, R. K. (1974). Predictability, constancy, and contingency 
#'             of periodic phenomena. Ecology, 55(5), 1148-1153.
calculate_colwell_indices <- function(x) {
  
  # Set NA handling
  if (all(is.na(x))) {
    return(c(Pc = NA, Pm = NA))
  }
  
  # --- 1. Data Structuring and Normalisation ---
  
  # The input time series (x) is 294 months long (24 full years + 6 months).
  # We need a matrix of 12 columns (Months) where each row is a Year.
  
  M_TOTAL <- length(x) # 294 months
  N_MONTHS <- 12       # Number of months in a year
  
  # Pad the data to create a full 25th row (25 * 12 = 300 months). 
  # The last 6 months of 2025 are already missing/NA, so we pad the last 6 months 
  # with NA to make the final matrix square (300 cells).
  x_padded <- c(x, rep(NA, N_MONTHS * ceiling(M_TOTAL / N_MONTHS) - M_TOTAL))
  
  # Reshape into a matrix: Rows = Years (25), Cols = Months (12)
  M_ym <- matrix(x_padded, 
                 ncol = N_MONTHS, 
                 byrow = TRUE)
  
  # --------------------------------------------------------------------------
  # 2. Proportion Calculation (Normalisation)
  # All calculations are based on proportions relative to the total sum.
  # --------------------------------------------------------------------------
  
  # Total sum of the observation across all time (used for normalization)
  # We use the original vector x here to avoid summing the padded NAs.
  X_total <- sum(x, na.rm = TRUE)
  
  # If the total observation is near zero (e.g., permanent desert), return NA
  if (X_total < 1e-6) {
    return(c(Pc = NA, Pm = NA))
  }
  
  # P_ij: Proportion of observation at time i (month) and year j.
  # This matrix is the normalized version of M_ym.
  P_ij <- M_ym / X_total
  
  # P_i: Proportion of observation in month i (Mean proportion for each month across all years). 
  # This is the mean of each column.
  P_i <- colMeans(P_ij, na.rm = TRUE)
  
  # P_j: Proportion of observation in year j (Total proportion for each year).
  # This is the sum of each row.
  P_j <- rowSums(P_ij, na.rm = TRUE)
  
  
  # --------------------------------------------------------------------------
  # 3. Apply Colwell's Formulae (using Entropy/Log-based approach)
  # The entropy measure H = -sum(p * log(p)).
  # --------------------------------------------------------------------------
  
  # Define helper function for entropy (H) calculation
  calculate_H <- function(probs) {
    # Exclude non-positive values as log(0) is undefined and log(NA) is NA.
    probs <- probs[probs > 0 & !is.na(probs)]
    if (length(probs) == 0) return(0)
    # H = -sum(p * log(p))
    return(-sum(probs * log(probs)))
  }
  
  # Maximum entropy (H_max) for 12 months (log(12))
  H_max <- log(N_MONTHS) 
  
  # --- a. Total Predictability (P) ---
  # P = 1 - (H_total / H_max)
  # H_total is the entropy of the proportions P_i (mean month proportions).
  H_total <- calculate_H(P_i)
  
  P <- 1 - (H_total / H_max)
  
  # --- b. Constancy (Pc) ---
  # Pc = 1 - (H_month / H_max)
  # H_month measures the average monthly variation regardless of the specific month.
  
  # H_month is the mean entropy of the P_j vector (variation across years for each month).
  H_month_vector <- apply(P_ij, MARGIN = 2, FUN = calculate_H)
  H_month <- mean(H_month_vector, na.rm = TRUE)
  
  # This calculation is slightly non-standard due to the partial year, 
  # but follows the original decomposition logic:
  Pc <- 1 - (H_month / H_max)
  
  # --- c. Contingency (Pm) ---
  # Pm is the portion of P that is not explained by Pc.
  Pm <- P - Pc
  
  # --------------------------------------------------------------------------
  # 4. Final Validation and Return
  # --------------------------------------------------------------------------
  
  # Colwell's indices range from 0 to 1. Minor numerical noise might put them 
  # slightly outside, so we constrain the values.
  Pc <- max(0, min(1, Pc))
  Pm <- max(0, min(1, Pm))
  
  return(c(Pc = Pc, Pm = Pm))
}