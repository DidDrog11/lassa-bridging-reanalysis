library(terra)
library(tidyterra)
library(ggplot2)
library(ggfx)
library(rnaturalearth)
library(maps)
library(readr)
library(dplyr)
library(tidyr)

# Load spatial assets
bbox <- ext(-18, 16, 4, 26)

admin_clean <- ne_countries(scale = "medium", returnclass = "sv") |>
  crop(bbox)
pts <- read_rds(here::here("data", "processed", "rodent_data_mapped.rds")) |>
  vect(geom = c("lon", "lat"), crs = crs(admin_clean)) |>
  filter(observed >= 1)

# Get outer boundary line string of the entire region
coast_line <- as.lines(aggregate(admin_clean))
coast_box <- ext(-18, 16, 4, 15)
clean_coast <- crop(coast_line, coast_box)
void_buf <- crop(buffer(clean_coast, width = 75000), admin_clean)

# Generate a 75km buffer and mask it strictly to the landmass
void_buf <- crop(buffer(clean_coast, width = 75000), admin)

cities_filtered <- maps::world.cities |> 
  subset(long >= -18 & long <= 16 & lat >= 4 & lat <= 26 & pop > 250000) |> 
  group_by(country.etc) |> 
  arrange(desc(pop)) |> 
  filter(!(country.etc == "Nigeria" & row_number() > 3)) |> 
  ungroup()

large_cities <- vect(cities_filtered, geom = c("long", "lat"), crs = "EPSG:4326")

# Simple dark theme layout for presentation slide consistency
theme_dark_slide <- function() {
  theme_void() + 
    theme(
      plot.background = element_rect(fill = "#05080f", colour = NA), 
      panel.background = element_rect(fill = "#05080f", colour = NA),
      legend.position = "inside",
      legend.position.inside = c(0.82, 0.75), 
      legend.text = element_text(colour = "#ffffff", size = 11),
      legend.title = element_blank(),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.spacing.y = unit(0.2, "cm")
    )
}

# Construct the neon density map
p1 <- ggplot() +
  geom_spatvector(data = admin_clean, fill = "#0a0f1d", colour = "#161f30", linewidth = 1.2) +
  with_outer_glow(geom_spatvector(data = void_buf, fill = "#ff0033", alpha = 0.15, colour = NA), colour = "#ff0033", sigma = 25, expand = 8) +
  with_outer_glow(geom_spatvector(data = pts, aes(colour = "Rodent Sampling Locations"), size = 0.15, alpha = 0.35), colour = "#00e5ff", sigma = 3) +
  with_outer_glow(geom_spatvector(data = large_cities, aes(colour = "Large Cities (>250k)"), size = 0.8, alpha = 0.8), colour = "#fff200", sigma = 2) +
  geom_spatvector_text(data = large_cities, aes(label = name), colour = "#ffffff", size = 2.8, fontface = "bold", check_overlap = TRUE, nudge_y = 0.35) +
  scale_colour_manual(values = c("Rodent Sampling Locations" = "#00e5ff", "Large Cities (>250k)" = "#fff200")) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1, shape = 16))) +  annotate("text", x = 15, y = 25, label = "SPATIAL BIAS IN WEST AFRICAN RODENT SAMPLING",
           colour = "#ffffff", size = 5.5, fontface = "bold", hjust = 1, vjust = 1) +
  annotate("text", x = 15, y = 23.8, label = "The Urban 'Sampling Void' vs. Neon Population Centres",
           colour = "#627791", size = 4.5, fontface = "italic", hjust = 1, vjust = 1) +
  theme_dark_slide()

ggsave(here::here("results", "presentation", "presentation_void_map.png"),
       plot = p1, width = 11, height = 6, dpi = 300)

# Load the previously calculated and normalised radial profiles
normalised_df <- read_rds(here::here("data", "processed", "normalised_profile_data.rds"))

# Dedicated dark theme for plots with axes
theme_dark_profile <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "#05080f", colour = NA),
      panel.background = element_rect(fill = "#05080f", colour = NA),
      panel.grid.major = element_line(colour = "#1a2333", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      text = element_text(colour = "#ffffff"),
      axis.text = element_text(colour = "#8b9eb3"),
      axis.title = element_text(colour = "#ffffff", face = "bold", margin = margin(t = 10, r = 10)),
      plot.title = element_text(colour = "#ffffff", face = "bold", size = 18),
      plot.subtitle = element_text(colour = "#627791", face = "italic", margin = margin(b = 15)),
      strip.text = element_text(colour = "#ffffff", face = "bold", size = 14),
      legend.position = "bottom",
      legend.text = element_text(colour = "#ffffff"),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )
}

plot_urban_paradox_dark <- function(data, title_text) {
  peak_stats <- data |> 
    filter(Variable == "Incidence") |> 
    group_by(City) |> 
    arrange(desc(Mean)) |> 
    slice_head(n = 1) |> 
    summarise(
      Peak_Dist = Distance,
      Peak_Val_Raw = round(Mean, 1),
      Peak_Norm_Y = Norm_Mean
    ) |> 
    mutate(
      Text_X = ifelse(Peak_Dist < 20, 35, 10),
      Text_Y = 0.85
    )
  
  ggplot(data, aes(x = Distance, y = Norm_Mean, colour = Variable_Label)) +
    geom_point(data = data |> filter(Variable == "Incidence"), 
               aes(x = Distance, y = Norm_Mean), 
               colour = "#fff200", alpha = 0.15, size = 1) +
    with_outer_glow(
      geom_smooth(aes(linetype = Variable_Label), linewidth = 1.2, se = FALSE, span = 0.3),
      colour = "white", sigma = 2, expand = 1
    ) +
    geom_segment(data = peak_stats,
                 aes(x = Text_X, xend = Peak_Dist, y = Text_Y, yend = Peak_Norm_Y),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                 colour = "#fff200", linewidth = 0.6, alpha = 0.9, inherit.aes = FALSE) +
    geom_text(data = peak_stats,
              aes(x = Text_X, y = Text_Y + 0.05, 
                  label = paste0("Peak: ", Peak_Val_Raw, "\n(", round(Peak_Dist,0), "km)")),
              colour = "#fff200", fontface = "bold", size = 3.5, lineheight = 0.9, 
              inherit.aes = FALSE) +
    scale_colour_manual(values = c("Ecological Risk (Dx)" = "#ff5500",  
                                   "Urban Shield (NTL)"   = "#00e5ff",  
                                   "Predicted Incidence"  = "#fff200")) + 
    scale_linetype_manual(values = c("Ecological Risk (Dx)" = "solid", 
                                     "Urban Shield (NTL)"   = "dotted",
                                     "Predicted Incidence"  = "solid")) +
    facet_wrap(~City, ncol = 3) +
    labs(title = title_text,
         subtitle = "Smoothed radial profiles showing the decoupling of risk and incidence.",
         x = "Distance from City Centre (km)",
         y = "Normalised Intensity",
         colour = NULL, linetype = NULL) +
    theme_dark_profile()
}

main_cities <- c("Lagos", "Tamale", "Jos")
plot_data_main <- normalised_df |> 
  filter(City %in% main_cities) |>
  mutate(City = factor(City, levels = main_cities))

p_main_dark <- plot_urban_paradox_dark(plot_data_main, "SPATIAL DECOUPLING: HAZARD VS. INCIDENCE")

ggsave(here::here("results", "presentation", "presentation_urban_paradox.png"), 
       plot = p_main_dark, width = 12, height = 5.5, dpi = 300, bg = "#05080f")

library(patchwork)

# ==============================================================================
# PART 3: TRUE-GEOGRAPHY NEON CITY ZOOMS
# ==============================================================================

# Ensure the original raster stack is loaded (from your earlier script)
Dm <- terra::rast(here(maps_dir, "Dm_Mastomys_natalensis.tif"))
Dl <- terra::rast(here(maps_dir, "Dl_Pathogen_Layer_Ensemble.tif"))
Dx <- Dm * Dl
names(Dx) <- "Risk_Dx"


# Shield (NTL)
stack_path <- here(proc_dir, "final_predictor_stack_scaled.tif")
NTL <- terra::rast(stack_path)[["NTL"]]
names(NTL) <- "Shield_NTL"

# Outcome (Incidence)
# Loading the final masked/calculated incidence
inc_path <- here(maps_dir, "Incidence_Full_Custom_Lambda0.03.tif")
if(file.exists(inc_path)) {
  Incidence <- terra::rast(inc_path)
} else { stop("Incidence raster not found.") }
names(Incidence) <- "Incidence"

# Population (for context if needed, strictly we use Incidence which includes Pop)
# But plotting Pop density helps explain the peri-urban peak
pop_path <- here(proc_dir, "pop_count_2020_05deg.tif")
Pop <- terra::rast(pop_path)
names(Pop) <- "Population"

# Stack all for extraction
# Ensure alignment (resample NTL/Pop to Dx if needed, usually they match)
if(!compareGeom(Dx, NTL, stopOnError=FALSE)) NTL <- resample(NTL, Dx)
if(!compareGeom(Dx, Incidence, stopOnError=FALSE)) Incidence <- resample(Incidence, Dx)
if(!compareGeom(Dx, Pop, stopOnError=FALSE)) Pop <- resample(Pop, Dx)

City_Stack <- c(Dx, NTL, Incidence, Pop)

plot_real_geography_dark <- function(city_name, lon, lat, raster_stack) {
  
  # 1. Define a 60km bounding box around the city center
  city_pt <- vect(matrix(c(lon, lat), ncol=2), crs="EPSG:4326")
  city_box <- buffer(city_pt, width = 60000) |> ext()
  
  # 2. Fetch HIGH-RESOLUTION coastlines and crop
  local_admin <- ne_countries(scale = "large", returnclass = "sf") |> 
    vect() |> 
    crop(city_box)
  
  # 3. Crop the raw rasters
  local_stack <- crop(raster_stack, city_box)
  
  # 4. SMOOTHING: Create a high-res template (e.g., 0.002 degrees)
  hires_template <- rast(city_box, res = 0.002, crs = crs(local_stack))
  smooth_stack <- resample(local_stack, hires_template, method = "bilinear")
  
  # 5. Mask to landmass 
  smooth_stack <- mask(smooth_stack, local_admin)
  
  # 6. DISTANCE CONTOURS: Calculate radar rings
  dist_rast <- distance(hires_template, city_pt)
  dist_contours <- as.contour(dist_rast, levels = seq(10000, 50000, by = 10000))
  
  # 7. Normalise the Shield locally for the Alpha channel
  scale_raster <- function(r) { 
    rng <- minmax(r)
    if(is.na(rng[1]) || rng[1] == rng[2]) return(r) 
    (r - rng[1]) / (rng[2] - rng[1]) 
  }
  
  ntl_norm <- scale_raster(smooth_stack[["Shield_NTL"]])
  dx_norm  <- scale_raster(smooth_stack[["Risk_Dx"]])
  
  # 8. Extract Data Frames for ggplot
  ntl_df <- as.data.frame(ntl_norm, xy = TRUE, na.rm = TRUE)
  dx_df  <- as.data.frame(dx_norm, xy = TRUE, na.rm = TRUE)
  
  # 9. INCIDENCE FIX: Calculate contours and split by level
  inc_vals <- values(smooth_stack[["Incidence"]], mat = FALSE, na.rm = TRUE)
  inc_breaks <- quantile(inc_vals[inc_vals > 0], probs = c(0.85, 0.95), na.rm = TRUE)
  
  if(length(unique(inc_breaks)) < 2) {
    inc_breaks <- c(max(inc_vals, na.rm=T) * 0.85, max(inc_vals, na.rm=T) * 0.95)
  }
  
  inc_contour <- as.contour(smooth_stack[["Incidence"]], levels = inc_breaks)
  
  # Split contours by level for distinct styling
  levels_present <- sort(unique(inc_contour$level))
  inc_low  <- inc_contour[inc_contour$level == levels_present[1], ]
  
  # Safety catch: Ensure the 95th percentile contour actually generated a geometry
  has_high <- length(levels_present) > 1
  if(has_high) {
    inc_high <- inc_contour[inc_contour$level == levels_present[2], ]
  }
  
  # 10. Plotting
  p <- ggplot() +
    geom_spatvector(data = local_admin, fill = "#05080f", colour = "#1a2333", linewidth = 0.5) +
    
    geom_raster(data = dx_df, aes(x = x, y = y, fill = Risk_Dx)) +
    scale_fill_gradientn(
      colours = c("#1a0b00", "#8c2d04", "#ff5500"), 
      na.value = "transparent", 
      guide = "none"
    ) +
    
    geom_raster(data = ntl_df, aes(x = x, y = y, alpha = Shield_NTL), fill = "#00e5ff") +
    scale_alpha_continuous(range = c(0, 0.85), na.value = 0, guide = "none") +
    
    geom_spatvector(data = dist_contours, colour = "#4a5a75", linewidth = 0.4, linetype = "dashed", alpha = 0.6) +
    
    # Lower Incidence (85th percentile) - Glowing Yellow
    with_outer_glow(
      geom_spatvector(data = inc_low, colour = "#fff200", linewidth = 1),
      colour = "#fff200", sigma = 3, expand = 1
    ) +
    
    # Higher Incidence (95th percentile) - Glowing Neon Red
    (if(has_high) {
      with_outer_glow(
        geom_spatvector(data = inc_high, colour = "#ff0033", linewidth = 1.5),
        colour = "#ff0033", sigma = 5, expand = 2
      )
    } else {
      NULL
    }) +
    
    geom_point(aes(x = lon, y = lat), colour = "white", shape = 3, size = 3) +
    
    labs(title = toupper(city_name)) +
    theme_dark_profile(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
  return(p)
}

# Generate the plots (Ensure City_Stack is in your environment)
p_lagos  <- plot_real_geography_dark("Lagos", 3.3792, 6.5244, City_Stack)
p_tamale <- plot_real_geography_dark("Tamale", -0.8393, 9.4008, City_Stack)
p_jos    <- plot_real_geography_dark("Jos", 8.8853, 9.8965, City_Stack)

# Stitch together
combined_geography <- p_lagos | p_tamale | p_jos

final_geography_plot <- combined_geography + 
  plot_annotation(
    title = "ANISOTROPIC SPATIAL DECOUPLING",
    subtitle = "Geography shapes the peri-urban transmission fringe. (Cyan = Shield, Orange = Hazard, Yellow = Incidence)",
    theme = theme(
      plot.background = element_rect(fill = "#05080f", colour = NA),
      plot.title = element_text(colour = "#ffffff", face = "bold", size = 20),
      plot.subtitle = element_text(colour = "#627791", face = "italic", size = 14)
    )
  )

ggsave(here::here("results", "presentation", "presentation_geography_zooms.png"), 
       plot = final_geography_plot, width = 14, height = 6, dpi = 300, bg = "#05080f")
