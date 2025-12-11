################################################################################
## SCRIPT: 2_host_layer_interpretation.R
## PURPOSE: Compare the IMSOM Host Layer (Custom) vs Basinski (2021).
##          Demonstrate the "Infilling" of urban reservoir niches.
## INPUTS:  Dm_Mastomys_natalensis.tif, basinski_host_raster.tif, 
##          pop_count_2020_05deg.tif
## OUTPUTS: Figures 2 and 3 and supplementary
################################################################################

source(here::here("packages.R"))

# Directories
maps_dir <- here("results", "maps")
figs_dir <- here("results", "figures")
model_dir <- here("results", "models")
data_dir <- here("data", "raw")
proc_dir <- here("data", "processed")

# 1. Load Host Rasters ----------------------------------------------------
# A. Custom Model (IMSOM)
Dm_custom <- terra::rast(here(maps_dir, "Dm_Mastomys_natalensis.tif"))
names(Dm_custom) <- "Custom"

# B. Basinski Model (Reference)
# Using the comparison raster generated/loaded in previous steps
bas_path <- here(data_dir, "basinski_comparison", "host_raster.tif")
if(file.exists(bas_path)) {
  Dm_basinski <- terra::rast(bas_path)
  # Resample to match Custom if needed
  if(!terra::compareGeom(Dm_custom, Dm_basinski, stopOnError=FALSE)) {
    message("Resampling Basinski raster to match grid...")
    Dm_basinski <- terra::resample(Dm_basinski, Dm_custom, method="bilinear")
  }
} else {
  stop("Basinski Host Raster not found. Cannot perform comparison.")
}
names(Dm_basinski) <- "Basinski"

iucn_path <- here(data_dir, "iucn", "data_0.shp")
if(file.exists(iucn_path)) {
  iucn_vect <- vect(iucn_path)
  # Create a binary mask raster (1 Inside, NA Outside)
  iucn_mask <- terra::rasterize(iucn_vect, Dm_custom, field=1, background=NA)
  
  # Apply Mask to BOTH models (Standardizing the comparison universe)
  Dm_custom   <- terra::mask(Dm_custom, iucn_mask)
  Dm_basinski <- terra::mask(Dm_basinski, iucn_mask)
} else { warning("IUCN Shapefile not found. Skipping mask.") }

# C. Human Population 
Pop_rast <- terra::rast(here(proc_dir, "pop_count_2020_05deg.tif"))
# Log-transform for better plotting (log10(Pop + 1))
LogPop <- log10(Pop_rast + 1)
names(LogPop) <- "Log_Pop"

# D. Context (West Africa Borders)
west_africa_countries <- countrycode::codelist |>
  filter(str_detect(region23, "Western Africa")) |>
  select(country = country.name.en, iso3c, iso2c) |>
  filter(!str_detect(iso3c, "CPV|SHN"))
west_africa_countries$country[west_africa_countries$country == "Côte d’Ivoire"] <- "Ivory Coast"

wa_vect <- ne_countries(scale = 10, country = west_africa_countries$country, returnclass = "sv")

# Target cities
standardised_cities <- read_rds(here(proc_dir, "standardised_cities.rds"))

# IMSOM model
imsom_path <- here(model_dir, "imsom_model_fit.rds")
if(file.exists(imsom_path)) {
  out_imsom <- readRDS(imsom_path)
} else { stop("IMSOM Model object not found.") }

# 2. Global Difference Map ------------------------------------------------
# Panel A: The IMSOM Host Map
p_host_map <- ggplot() +
  geom_spatraster(data = Dm_custom) +
  geom_spatvector(data = wa_vect, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", name = "P(Occ)", na.value = "transparent", limits = c(0, 1)) +
  theme_minimal() +
  labs(title = "A. Predicted Reservoir Niche (Biotic IMSOM)",
       x = NULL, y = NULL) +
  theme(legend.position = "right", legend.key.height = unit(1.2, "cm"))

# Panel B: The Difference Map
Dm_diff <- Dm_custom - Dm_basinski
names(Dm_diff) <- "Difference"

p_diff <- ggplot() +
  geom_spatraster(data = Dm_diff) +
  geom_spatvector(data = wa_vect, fill = NA, colour = "black", linewidth = 0.2) +
  # Purple/Green Divergent Palette
  scale_fill_gradient2(low = "#1B7837", mid = "white", high = "#762A83", midpoint = 0,
                       name = expression(Delta ~ P), na.value = "transparent") +
  theme_minimal() +
  labs(title = "B. Model Divergence",
       subtitle = "Purple: Biotic IMSOM Higher | Green: Climatic SDM Higher",
       x = NULL, y = NULL) +
  theme(legend.position = "right", legend.key.height = unit(1.2, "cm"))

# Combine and Save Figure 2
fig2 <- p_host_map / p_diff

ggsave(here(figs_dir, "Fig2_Host_Maps.png"), fig2, width = 8, height = 10, bg="white")

# 3. Urban Transects -----------------------------
buffer_size <- 0.6 # Degrees (~60km box)
plot_list <- list()

for(i in 1:nrow(standardised_cities)) {
  city <- standardised_cities[i,]
  # Crop Rasters
  e <- ext(city$lon - buffer_size, city$lon + buffer_size, 
           city$lat - buffer_size, city$lat + buffer_size)
  r_cust <- terra::crop(Dm_custom, e)
  r_base <- terra::crop(Dm_basinski, e)
  # Create Distance Raster for Contours (Distance from City Center in km)
  dist_rast <- r_cust
  # Set values to coordinate distance
  xy <- terra::xyFromCell(dist_rast, 1:ncell(dist_rast))
  dists <- sqrt((xy[,1] - city$lon)^2 + (xy[,2] - city$lat)^2) * 111
  values(dist_rast) <- dists
  names(dist_rast) <- "Dist_km"
  
  # Convert to DF
  df_cust <- as.data.frame(r_cust, xy = TRUE, na.rm = FALSE) |> 
    mutate(Model = "Biotic IMSOM") |> 
    rename(Val = Custom)
  df_base <- as.data.frame(r_base, xy = TRUE, na.rm = FALSE) |> 
    mutate(Model = "Climate SDM") |> 
    rename(Val = Basinski)
  df_city <- bind_rows(df_base, df_cust)
  
  df_city <- df_city |>
    mutate(Dist_km = sqrt((x - city$lon)^2 + (y - city$lat)^2) * 111)
  
  # Plot
  p <- ggplot(df_city, aes(x = x, y = y)) +
    geom_raster(aes(fill = Val)) +
    geom_contour(aes(z = Dist_km), breaks = c(20, 40), 
                 color = "cyan", linetype = "dotted", linewidth = 0.6, alpha = 0.8) +
    facet_wrap(~Model) +
    scale_fill_viridis_c(option = "magma", limits = c(0,1), name = "P(Occ)", 
                         na.value = "#D1E5F0") +
    annotate("point", x = city$lon, y = city$lat, colour = "cyan", shape = 3, size = 3, stroke = 1.5) +
    annotate("segment", 
             x = e[1] + 0.8*buffer_size, xend = e[1] + 0.8*buffer_size + 0.18, 
             y = e[3] + 0.1*buffer_size, yend = e[3] + 0.1*buffer_size, 
             color = "cyan", linewidth = 1) +
    annotate("text", 
             x = e[1] + 0.8*buffer_size + 0.09, 
             y = e[3] + 0.1*buffer_size + 0.05, 
             label = "20km", color = "cyan", size = 2.5, fontface = "bold") +
    coord_fixed() +
    theme_void() +
    labs(title = paste0(city$name, " (", city$type, ")")) +
    theme(legend.position = "none", 
          strip.text = element_text(face = "bold", size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plot_list[[city$name]] <- p
}

p_zooms_main <- (plot_list[["Lagos"]] / plot_list[["Tamale"]] / plot_list[["Jos"]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# All selected cities to supplementary
p_supp_1 <- wrap_plots(plot_list[1:6], ncol = 2)
p_supp_2 <- wrap_plots(plot_list[7:12], ncol = 2)

ggsave(here(figs_dir, "Supp_Urban_Zooms_Host_1.png"), p_supp_1, width = 10, height = 12, bg="white")
ggsave(here(figs_dir, "Supp_Urban_Zooms_Host_2.png"), p_supp_2, width = 10, height = 12, bg="white")

# 4. Response Curve: Occupancy vs Population Density ----------------------
# Load Uncertainty Raster
Dm_sd <- terra::rast(here(maps_dir, "Dm_Mastomys_natalensis_SD.tif"))
names(Dm_sd) <- "Custom_SD"

stack_response <- c(Dm_custom, Dm_sd, Dm_basinski, LogPop)
set.seed(123)
samp <- terra::spatSample(stack_response, size = 100000, method = "regular", na.rm = TRUE)
samp_df <- as.data.frame(samp)

# Fit GAMs for Smooth Trends
# Fit Trends for Mean Occupancy
gam_custom   <- gam(Custom ~ s(Log_Pop, bs = "cs"), data = samp_df)
gam_basinski <- gam(Basinski ~ s(Log_Pop, bs = "cs"), data = samp_df)

# Fit Trend for Uncertainty
gam_sd_custom <- gam(Custom_SD ~ s(Log_Pop, bs = "cs"), data = samp_df)

# Generate Prediction Grid
# Create a smooth sequence from min to max population
pred_grid <- data.frame(Log_Pop = seq(min(samp_df$Log_Pop), max(samp_df$Log_Pop), length.out = 200))

# Predict Means
pred_grid$Custom_Mean   <- predict(gam_custom, newdata = pred_grid, type = "response")
pred_grid$Basinski_Mean <- predict(gam_basinski, newdata = pred_grid, type = "response")

# Predict Uncertainty (SD)
pred_grid$Custom_SD_Pred <- predict(gam_sd_custom, newdata = pred_grid, type = "response")

# Reshape for Plotting
plot_df <- bind_rows(
  pred_grid |> 
    select(Log_Pop, Occupancy = Custom_Mean) |> 
    mutate(Model = "Biotic IMSOM", 
           # Attach the smoothed SD for the ribbon
           SD = pred_grid$Custom_SD_Pred),
  pred_grid |> 
    select(Log_Pop, Occupancy = Basinski_Mean) |> 
    mutate(Model = "Climate SDM", 
           SD = NA) # Basinski has no SD raster
)

# Plot
p_response <- ggplot(plot_df, aes(x = Log_Pop, y = Occupancy, colour = Model, fill = Model)) +
  geom_ribbon(data = filter(plot_df, Model == "Biotic IMSOM"),
              aes(ymin = pmax(0, Occupancy - SD), 
                  ymax = pmin(1, Occupancy + SD)), 
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("Climate SDM" = "#0072B2", "Biotic IMSOM" = "#D55E00")) +
  scale_fill_manual(values = c("Climate SDM" = "#0072B2", "Biotic IMSOM" = "#D55E00")) +
  theme_minimal() +
  labs(title = "Host Tolerance to Urbanisation",
       x = "Log10 Human Population Density",
       y = "Predicted Occupancy Probability") +
  theme(legend.position = c(0.85, 0.1),
        legend.background = element_rect(fill = "white", colour = NA))

# 5. Save -----------------------------------------------------------------
# Map + Response
fig3 <- cowplot::plot_grid(p_zooms_main, p_response,
                           labels = c('A.', 'B.'),
                           label_size = 12,
                           label_fontface = "plain",
                           rel_widths = c(1, 1.4))
ggsave(here(figs_dir, "Fig3_Host_Mechanisms.png"), fig3, width = 12, height = 8, bg="white")

# 6. IMSOM coefficients ---------------------------------------------------
# Extract Beta Samples
beta_mat <- as.matrix(out_imsom$beta.samples)

# Summarise
beta_summary <- data.frame(
  Parameter = colnames(beta_mat),
  Mean = apply(beta_mat, 2, mean),
  Lower = apply(beta_mat, 2, quantile, probs = 0.025),
  Upper = apply(beta_mat, 2, quantile, probs = 0.975)
) |>
  tibble::rownames_to_column("row_id") |>
  filter(!grepl("\\(Intercept\\)", Parameter)) |>
  separate(Parameter, into = c("Covariate", "Species"), sep = "-", extra = "merge")

# Formatting for Plot
covariate_labels <- c(
  "Tmu"           = "Mean Temperature",
  "Pmu"           = "Mean Precipitation",
  "Nmu"           = "Mean Greenness (NDVI)",
  "Elev"          = "Elevation",
  "Pop"           = "Human Population",
  "I(Pop^2)"      = "Human Population^2 (Quadratic)",
  "LC_12_Density" = "Cropland Density"
)

plot_data <- beta_summary |>
  mutate(Species = gsub("_", " ", Species),
         Highlight = case_when(
           Species == "Mastomys natalensis" ~ "Target",
           Species %in% c("Rattus rattus", "Mus musculus") ~ "Invasive",
           TRUE ~ "Native"
         )) |>
  mutate(Covariate_Clean = recode(Covariate, !!!covariate_labels)) |>
  mutate(Covariate_Clean = coalesce(Covariate_Clean, Covariate))

# Generate Plot
p_coef <- ggplot(plot_data, aes(x = Mean, y = reorder(Species, Mean), colour = Highlight)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~Covariate_Clean, scales = "free_x", nrow = 2, labeller = label_wrap_gen(width = 25)) +
  scale_colour_manual(values = c("Target" = "#D55E00", "Invasive" = "#0072B2", "Native" = "grey60")) +
  theme_bw(base_size = 14) +
  labs(title = "Predictors of Rodent Occupancy",
       x = "Effect Size (Logit Scale)", 
       y = NULL, 
       colour = "Group") +
  theme(strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        axis.text.y = element_text(face = "italic"))

# Save
ggsave(here(figs_dir, "Supp_Host_Predictors.png"), p_coef, width = 14, height = 10, bg="white")
