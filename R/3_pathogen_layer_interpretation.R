################################################################################
## SCRIPT: 3_pathogen_layer_interpretation.R
## PURPOSE: Analyze the drivers of Lassa Virus Prevalence (Dl).
##          Test the "Biotic Amplification" hypothesis (Rattus effect).
## INPUTS:  pathogen_model_fit_full.rds (BRT Model), 
##          Dl_Pathogen_Layer_Ensemble.tif (Raster)
## OUTPUTS: Figures 4 and 5
################################################################################

source(here::here("packages.R"))

# Directories
maps_dir <- here("results", "maps")
figs_dir <- here("results", "figures")
proc_dir <- here("data", "processed")
model_dir <- here("results", "models")
data_dir <- here("data", "raw")

if(!dir.exists(figs_dir)) dir.create(figs_dir)

# 1. Load BRT Model and Raster --------------------------------------------
# A. BRT Model
brt_path <- here(model_dir, "pathogen_model_fit_full.rds")
if(file.exists(brt_path)) brt_model <- readRDS(brt_path) else stop("BRT Model not found.")

# B. Custom Pathogen Raster (Dl)
Dl_custom <- terra::rast(here(maps_dir, "Dl_Pathogen_Layer_Ensemble.tif"))
names(Dl_custom) <- "Custom"

# C. Basinski Pathogen Raster (Reference)
bas_dl_path <- here(data_dir, "basinski_comparison", "pathogen_raster.tif")
if(file.exists(bas_dl_path)) {
  Dl_basinski <- terra::rast(bas_dl_path)
  # Align/Resample
  if(!terra::compareGeom(Dl_custom, Dl_basinski, stopOnError=FALSE)) {
    Dl_basinski <- terra::resample(Dl_basinski, Dl_custom, method="bilinear")
  }
} else { stop("Basinski Pathogen raster not found.") }
names(Dl_basinski) <- "Basinski"

# D. IUCN Mask (The Biological Filter)
iucn_path <- here(data_dir, "iucn", "data_0.shp")
if(file.exists(iucn_path)) {
  iucn_vect <- vect(iucn_path)
  iucn_mask <- terra::rasterize(iucn_vect, Dl_custom, field=1, background=NA)
  
  # Apply Mask to BOTH models
  Dl_custom   <- terra::mask(Dl_custom, iucn_mask)
  Dl_basinski <- terra::mask(Dl_basinski, iucn_mask)
} else { warning("IUCN Shapefile not found. Skipping mask.") }

# E. Context (Borders)
west_africa_countries <- countrycode::codelist |>
  filter(str_detect(region23, "Western Africa")) |>
  select(country = country.name.en, iso3c, iso2c) |>
  filter(!str_detect(iso3c, "CPV|SHN"))
west_africa_countries$country[west_africa_countries$country == "Côte d’Ivoire"] <- "Ivory Coast"

wa_vect <- ne_countries(scale = 10, country = west_africa_countries$country, returnclass = "sv")

Pop_rast <- terra::rast(here(proc_dir, "pop_count_2020_05deg.tif"))
# Log-transform for better plotting (log10(Pop + 1))
LogPop <- log10(Pop_rast + 1)
names(LogPop) <- "Log_Pop"

# 2. Variable Importance Plot ---------------------------------------------
clean_labels <- c(
  "Nmin"          = "Greenness Persistence (NDVI Min)",
  "Dm_Rattus"     = "Invasive Occupancy (R. rattus)",
  "Pmax"          = "Max Precipitation",
  "Tmu"           = "Mean Temperature",
  "Pc"            = "Precipitation Constancy",
  "Nmax"          = "Peak Greenness",
  "Pm"            = "Precipitation Contingency",
  "Ndur"          = "Greenness Duration",
  "LC_9_Density"  = "Savanna Density",
  "Dm_Mus"        = "Invasive Occupancy (M. musculus)",
  "Pmu"           = "Mean Precipitation",
  "Nmu"           = "Mean Greenness",
  "LC_10_Density" = "Grassland Density",
  "NTL"           = "Nighttime Lights (Urban)",
  "LC_12_Density" = "Cropland Density"
)

# Extract influence
summary_df <- summary(brt_model, plotit = FALSE) |>
  arrange(desc(rel.inf)) |>
  mutate(var_clean = recode(var, !!!clean_labels)) |>
  mutate(Highlight = case_when(
    var %in% c("Dm_Rattus", "Dm_Mus") ~ "Biotic (Invasive)",
    var == "NTL" ~ "Anthropogenic",
    TRUE ~ "Environmental"
  ))

p_varimp <- ggplot(summary_df, aes(x = reorder(var_clean, rel.inf), y = rel.inf, fill = Highlight)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Biotic (Invasive)" = "#D55E00", 
                               "Environmental" = "gray70",
                               "Anthropogenic" = "black")) +
  theme_minimal() +
  labs(title = "A. Predictors of LASV Prevalence",
       x = NULL, y = "Relative Influence (%)") +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank())

# 3. Partial Dependence Plot -------------------------------
predictors <- tribble(
  ~var, ~Category, ~Label,
  "Nmin", "Environmental", "Greenness Persistence (NDVI Min)",
  "Pmax", "Environmental", "Max Precipitation",
  "Tmu", "Environmental", "Mean Temperature",
  "Pc", "Environmental", "Precipitation Constancy",
  "Nmax", "Environmental", "Peak Greenness",
  "Pm", "Environmental", "Precipitation Contingency",
  "Ndur", "Environmental", "Greenness Duration",
  "LC_9_Density", "Environmental", "Savanna Density",
  "Pmu", "Environmental", "Mean Precipitation",
  "Nmu", "Environmental", "Mean Greenness",
  "LC_10_Density", "Environmental", "Grassland Density",
  "LC_12_Density", "Environmental", "Cropland Density",
  "Dm_Rattus", "Biotic", "Invasive Occupancy (R. rattus)",
  "Dm_Mus", "Biotic", "Invasive Occupancy (M. musculus)"
)

get_pdp_df <- function(model, var_name, category, label) {
  pdp_raw <- plot.gbm(model, i.var = var_name, return.grid = TRUE)
  colnames(pdp_raw) <- c("x", "y")
  pdp_raw$Variable <- var_name
  pdp_raw$Category <- category
  pdp_raw$Label <- label
  return(pdp_raw)
}

all_pdps <- purrr::pmap_dfr(list(predictors$var, predictors$Category, predictors$Label), 
                            function(v, c, l) get_pdp_df(brt_model, v, c, l))
all_pdps$Label <- factor(all_pdps$Label, levels = summary_df$var_clean)

top_pdp_data <- all_pdps |> 
  filter(Variable %in% c("Nmin", "Dm_Rattus"))

p_nmin <- ggplot(top_pdp_data |> filter(Variable == "Nmin"), aes(x = x, y = y)) +
  geom_line(color = "#009E73", linewidth = 1.2) +
  labs(title = "B. Environmental Predictor",
       subtitle = "Greenness Persistence (NDVI Min)",
       x = "NDVI Minimum (Scaled)", 
       y = "Marginal Effect on Prevalence") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_rattus <- ggplot(top_pdp_data |> filter(Variable == "Dm_Rattus"), aes(x = x, y = y)) +
  geom_line(color = "#D55E00", linewidth = 1.2) +
  labs(title = "C. Biotic Predictor",
       subtitle = "Invasive Co-occurrence (R. rattus)",
       x = "R. rattus Occupancy Probability", 
       y = "Marginal Effect on Prevalence") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

fig4 <- p_varimp + (p_nmin / p_rattus) + 
  plot_layout(widths = c(1.2, 1))

ggsave(here(figs_dir, "Fig4_Pathogen_Predictors.png"), fig4, width = 10, height = 5, bg="white")

# 4. Interaction Plot: Rattus x Greenness (Nmin) --------------------------
pdp_2d <- gbm::plot.gbm(brt_model, i.var = c("Dm_Rattus", "Nmin"), return.grid = TRUE)

p_interact <- ggplot(pdp_2d, aes(x = Dm_Rattus, y = Nmin, fill = y)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "inferno", name = "Risk\n(Logit)") +
  theme_minimal() +
  labs(title = "Covariation of Predictors",
       x = "Invasive Occupancy (R. rattus)",
       y = "Greenness Persistence (NDVI Min)") +
  theme(legend.position = "right")

ggsave(here(figs_dir, "Supp_Interaction.png"), p_interact, width = 6, height = 5, bg="white")

# 5. The Hazard Map (Dl) --------------------------------------------------
# Panel A: Custom Prevalence Map (Masked)
p_prev_map <- ggplot() +
  geom_spatraster(data = Dl_custom) +
  geom_spatvector(data = wa_vect, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_viridis_c(option = "rocket", direction = -1, name = "Prevalence", na.value = "transparent", limits = c(0, 1)) +
  theme_minimal() +
  labs(title = "A. Predicted Pathogen Prevalence (Biotic-Informed)",
       subtitle = "Probability of LASV circulation given reservoir presence.",
       x = NULL, y = NULL) +
  theme(legend.position = "right", legend.key.height = unit(1.2, "cm"))

# Panel B: Difference Map
Dl_diff <- Dl_custom - Dl_basinski
names(Dl_diff) <- "Difference"

p_diff <- ggplot() +
  geom_spatraster(data = Dl_diff) +
  geom_spatvector(data = wa_vect, fill = NA, colour = "black", linewidth = 0.2) +
  scale_fill_gradient2(low = "#1B7837", mid = "white", high = "#762A83", midpoint = 0,
                       name = expression(Delta ~ P), na.value = "transparent") +
  theme_minimal() +
  labs(title = "B. Model Divergence",
       subtitle = "Purple: Biotic-Informed Higher | Green: Climatic Baseline Higher",
       x = NULL, y = NULL) +
  theme(legend.position = "right", legend.key.height = unit(1.2, "cm"))

# Assembly Fig 5
fig5 <- p_prev_map / p_diff
ggsave(here(figs_dir, "Fig5_Pathogen_Maps.png"), fig5, width = 8, height = 10, bg="white")

# 6. Mechanisms -----------------------------------------------------------
# Target cities
standardised_cities <- read_rds(here(proc_dir, "standardised_cities.rds"))

buffer_size <- 0.6 
plot_list <- list()

for(i in 1:nrow(standardised_cities)) {
  city <- standardised_cities[i,]
  # Crop Rasters
  e <- ext(city$lon - buffer_size, city$lon + buffer_size, 
           city$lat - buffer_size, city$lat + buffer_size)
  r_cust <- terra::crop(Dl_custom, e)
  r_base <- terra::crop(Dl_basinski, e)
  # Create Distance Raster for Contours (Distance from City Center in km)
  dist_rast <- r_cust
  # Set values to coordinate distance
  xy <- terra::xyFromCell(dist_rast, 1:ncell(dist_rast))
  dists <- sqrt((xy[,1] - city$lon)^2 + (xy[,2] - city$lat)^2) * 111
  values(dist_rast) <- dists
  names(dist_rast) <- "Dist_km"
  
  # Convert to DF
  df_cust <- as.data.frame(r_cust, xy = TRUE, na.rm = FALSE) |> 
    mutate(Model = "Biotic-Informed") |> 
    rename(Val = Custom)
  df_base <- as.data.frame(r_base, xy = TRUE, na.rm = FALSE) |> 
    mutate(Model = "Climatic Baseline") |> 
    rename(Val = Basinski)
  df_city <- bind_rows(df_base, df_cust)
  
  df_city <- df_city |>
    mutate(Dist_km = sqrt((x - city$lon)^2 + (y - city$lat)^2) * 111)
  
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

ggsave(here(figs_dir, "Supp_Urban_Zooms_Pathogen_1.png"), p_supp_1, width = 10, height = 12, bg="white")
ggsave(here(figs_dir, "Supp_Urban_Zooms_Pathogen_2.png"), p_supp_2, width = 10, height = 12, bg="white")

# Panel B: Response Curve (Prevalence vs Urbanisation)
dl_sd_path <- here(maps_dir, "Dl_Pathogen_Uncertainty.tif")
if(file.exists(dl_sd_path)) {
  Dl_sd <- terra::rast(dl_sd_path)
  if(exists("iucn_mask")) Dl_sd <- terra::mask(Dl_sd, iucn_mask)
} else {
  Dl_sd <- Dl_custom 
  values(Dl_sd) <- 0
}
names(Dl_sd) <- "Custom_SD"

set.seed(123)
samp <- terra::spatSample(c(Dl_custom, Dl_sd, Dl_basinski, LogPop), size = 100000, method = "regular", na.rm = TRUE)
samp_df <- as.data.frame(samp)

gam_custom   <- gam(Custom ~ s(Log_Pop, bs = "cs"), data = samp_df)
gam_basinski <- gam(Basinski ~ s(Log_Pop, bs = "cs"), data = samp_df)
gam_sd       <- gam(Custom_SD ~ s(Log_Pop, bs = "cs"), data = samp_df)

pred_grid <- data.frame(Log_Pop = seq(min(samp_df$Log_Pop), max(samp_df$Log_Pop), length.out = 200))
pred_grid$Custom <- predict(gam_custom, newdata=pred_grid, type="response")
pred_grid$Basinski <- predict(gam_basinski, newdata=pred_grid, type="response")
pred_grid$SD <- predict(gam_sd, newdata=pred_grid, type="response")

plot_df <- bind_rows(
  pred_grid |> select(Log_Pop, Prevalence=Custom, SD) |> mutate(Model="Biotic-Informed"),
  pred_grid |> select(Log_Pop, Prevalence=Basinski) |> mutate(Model="Climatic Baseline", SD=NA)
)

p_response <- ggplot(plot_df, aes(x = Log_Pop, y = Prevalence, colour = Model, fill = Model)) +
  geom_ribbon(data = filter(plot_df, Model == "Biotic-Informed"),
              aes(ymin = pmax(0, Prevalence - SD), ymax = pmin(1, Prevalence + SD)), 
              alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Climatic Baseline" = "#0072B2", "Biotic-Informed" = "#D55E00")) +
  scale_fill_manual(values = c("Climatic Baseline" = "#0072B2", "Biotic-Informed" = "#D55E00")) +
  theme_minimal() +
  labs(title = "Functional Response (Pathogen)",
       x = "Human Pop. Density (Log10)", y = "Predicted Prevalence") +
  theme(legend.position = c(0.8, 0.8), legend.background = element_rect(fill="white", colour=NA))

# Assembly Fig 6
fig6 <- cowplot::plot_grid(p_zooms_main, p_response,
                           labels = c('A.', 'B.'),
                           label_size = 12,
                           label_fontface = "plain",
                           rel_widths = c(1, 1.4))
ggsave(here(figs_dir, "Fig6_Pathogen_Mechanisms.png"), fig6, width = 12, height = 8, bg="white")
