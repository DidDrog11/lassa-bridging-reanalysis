# ---------------------------------------------------------------------------- #
# QC: Quality Control of Rodent Data
# ---------------------------------------------------------------------------- #

# Load necessary visualization packages
library(ggplot2)
library(tidyr)

# 1. Replicate Integrity Check -----------------------------------------------
# Goal: Ensure we have sites with >1 replicate (essential for estimating 'p')
# and that replicates aren't just 100 duplicates of the same day.

replicate_summary <- all_rodent_data_clean %>%
  group_by(cell_id) %>%
  summarise(
    n_replicates = n_distinct(replicate_id),
    n_species_detected = n_distinct(species[observed == 1]),
    total_surveys = n()
  )

message("--- Replicate Statistics (ArHa Data) ---")
print(summary(replicate_summary$n_replicates))

# Check: How many sites have repeated visits?
sites_with_repeats <- sum(replicate_summary$n_replicates > 1)
message(paste("Number of sites with > 1 temporal replicate:", sites_with_repeats))

if (sites_with_repeats < 10) {
  warning("Very few sites have temporal replicates! Detection probability (p) may be hard to estimate.")
}

# 2. Species Co-occurrence Matrix (Naive) ------------------------------------
# Goal: Before modeling residual covariance (rho), do we see raw overlap?
# We look at sites where species were *Observed* (1).

# Pivot to Wide Format (1 row per cell, columns for species presence)
site_sp_matrix <- all_rodent_data_clean %>%
  filter(observed == 1) %>%
  filter(species %in% target_species$rodent) %>% # Exclude Pseudo-Absences for this check
  distinct(cell_id, species) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = species, values_from = present, values_fill = 0) %>%
  select(-cell_id)

# Calculate simple co-occurrence count (Matrix Multiplication)
# 
co_occurrence_counts <- crossprod(as.matrix(site_sp_matrix))

message("--- Raw Species Co-occurrence Counts (Shared Sites) ---")
print(co_occurrence_counts)

# Quick check: Do M. natalensis and R. rattus overlap?
mn_rr_overlap <- co_occurrence_counts["Mastomys natalensis", "Rattus rattus"]
message(paste("Sites with BOTH M. natalensis and R. rattus:", mn_rr_overlap))


# 3. Spatial Distribution Plots (8 Species) ----------------------------------
# Goal: Visualize detection/non-detection for each target species.

# Get West Africa borders for context
wa_borders <- west_africa_ext |>
  aggregate("GID_0")

# Filter data for plotting
# We want to plot ArHa absences (0) and All Presences (1)
# We exclude the generic "Pseudo_Absence" species rows for the species-specific plots
plot_data <- all_rodent_data_clean %>%
  filter(species %in% target_species$rodent) %>%
  mutate(
    status = case_when(
      observed == 1 ~ "Present",
      observed == 0 ~ "Absent (ArHa)",
      TRUE ~ "Unknown"
    ),
    # Order factor for legend: Presence on top
    status = factor(status, levels = c("Present", "Absent (ArHa)"))
  )

# Generate the Faceted Plot
p_map <- ggplot() +
  # Base Map
  geom_sf(data = wa_borders, fill = "grey95", color = "grey80") +
  # Rodent Points
  geom_point(data = plot_data, 
             aes(x = lon, y = lat, color = status, shape = status),
             alpha = 0.7, size = 1.5) +
  # Colors: Red for Presence, Blue for Absence
  scale_color_manual(values = c("Present" = "#D55E00", "Absent (ArHa)" = "#0072B2")) +
  scale_shape_manual(values = c("Present" = 16, "Absent (ArHa)" = 4)) +
  # Facet by Species
  facet_wrap(~species, ncol = 4) +
  coord_sf(xlim = c(-17, 16), ylim = c(4, 20)) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.title = element_blank()
  ) +
  labs(title = "Spatial Distribution of Rodent Sampling by Species",
       subtitle = "Red = Detection, Blue = Non-detection (ArHa Surveys only)")

# Display Plot
print(p_map)

# Save Plot
ggsave(here(proc_dir, "misc", "qc_species_distributions.png"), p_map, width = 12, height = 6, bg = "white")