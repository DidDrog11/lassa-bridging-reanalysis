################################################################################
## SCRIPT: 5_urban_rural_incidence.R
## PURPOSE: Visualize the "Urban Paradox" via Radial Profiles.
##          Show divergence between Ecological Risk (Dx) and Incidence.
## INPUTS:  Dx (Risk), NTL (Shield), Incidence (Outcome), Pop (Density).
## OUTPUTS: Fig10_Urban_Paradox_Profiles.png, Supp_Fig2_All_City_Profiles.png
################################################################################

source(here::here("packages.R"))

# Directories
maps_dir <- here("results", "maps")
figs_dir <- here("results", "figures")
proc_dir <- here("data", "processed")
data_dir <- here("data", "raw")

if(!dir.exists(figs_dir)) dir.create(figs_dir)

# 1. Load Data ------------------------------------------------------------
# A. Rasters
# Risk (Dx) - The "Ecological Hazard"
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

# B. Define Cities (Main + Supp)
# Target cities
standardised_cities <- read_rds(here(proc_dir, "standardised_cities.rds"))

# 2. Extract Radial Profiles ----------------------------------------------
# Function to extract mean values in concentric rings
get_radial_profile <- function(city_name, lon, lat, raster_stack, max_dist_km = 50, step_km = 2) {
  
  pt <- vect(matrix(c(lon, lat), ncol=2), crs="EPSG:4326")
  dists <- seq(0, max_dist_km, by = step_km)
  profile_list <- list()
  
  for(i in 1:(length(dists)-1)) {
    
    # Geometry Logic
    buf_outer <- buffer(pt, width = dists[i+1] * 1000)
    
    if(dists[i] == 0) {
      ring <- buf_outer 
    } else {
      buf_inner <- buffer(pt, width = dists[i] * 1000)
      ring <- erase(buf_outer, buf_inner)
    }
    
    # Extract
    vals <- terra::extract(raster_stack, ring, ID=FALSE)
    
    if(nrow(vals) == 0) {
      # Return NA row if ring is empty (e.g. in ocean)
      res <- data.frame(Variable = names(raster_stack), Mean = NA, SD = NA)
    } else {
      # Summarise and Pivot
      res <- vals |>
        summarise(across(everything(), list(Mean = ~mean(., na.rm=TRUE), SD = ~sd(., na.rm=TRUE)))) |>
        pivot_longer(everything(), 
                     names_to = c("Variable", "Stat"), 
                     names_pattern = "^(.*)_(Mean|SD)$") |>
        pivot_wider(names_from = Stat, values_from = value)
    }
    
    profile_list[[i]] <- res |>
      mutate(City = city_name,
             Distance = (dists[i] + dists[i+1])/2)
  }
  bind_rows(profile_list)
}

# Loop through cities
all_profiles <- list()
for(i in 1:nrow(standardised_cities)) {
  cit <- standardised_cities[i,]
  all_profiles[[i]] <- get_radial_profile(cit$name, cit$lon, cit$lat, City_Stack)
}
profile_df <- bind_rows(all_profiles) |>
  left_join(standardised_cities[,c("name", "type")], by = c("City" = "name"))

# 3. Normalise Data (Scale 0-1) -------------------------------------------
# We normalise PER CITY to focus on the shape divergence
normalised_df <- profile_df |>
  group_by(City, Variable) |>
  mutate(Min_Val = min(Mean, na.rm=TRUE),
         Max_Val = max(Mean, na.rm=TRUE),
         Range = Max_Val - Min_Val,
         # Avoid division by zero if variable is constant
         Range = ifelse(Range == 0, 1, Range),
         Norm_Mean = (Mean - Min_Val) / Range,
         Norm_SD   = SD / Range,
         # Calculate Upper/Lower for Ribbons
         Lower = Norm_Mean - Norm_SD, 
         Upper = Norm_Mean + Norm_SD) |>
  ungroup() |>
  # Rename Variables for Plot
  mutate(Variable_Label = case_when(Variable == "Risk_Dx"    ~ "Ecological Risk (Dx)",
                                    Variable == "Shield_NTL" ~ "Urban Shield (NTL)",
                                    Variable == "Incidence"  ~ "Predicted Incidence"))


# 4. Plotting Function ----------------------------------------------------

plot_urban_paradox <- function(data, title_text) {
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
      Text_Y = 0.8 # Fixed height for alignment
    )
  
  ggplot(data, aes(x = Distance, y = Norm_Mean, colour = Variable_Label)) +
    geom_point(data = data |>
                 filter(Variable == "Incidence"), 
               aes(x = Distance, y = Norm_Mean, colour = Variable_Label), alpha = 0.3) +
    geom_smooth(aes(linetype = Variable_Label), size = 1.5, se = FALSE, span = 0.3) +
    geom_segment(data = peak_stats,
                 aes(x = Text_X, xend = Peak_Dist, 
                     y = Text_Y, yend = Peak_Norm_Y),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                 colour = "black", size = 0.6, alpha = 0.7, inherit.aes = FALSE) +
    geom_text(data = peak_stats,
              aes(x = Text_X, y = Text_Y, 
                  label = paste0("Peak: ", Peak_Val_Raw, "\n(", round(Peak_Dist,0), "km)")),
              colour = "black", fontface = "bold", size = 3.2, lineheight = 0.9, 
              inherit.aes = FALSE) +
    scale_colour_manual(values = c("Ecological Risk (Dx)" = "#D55E00", 
                                  "Urban Shield (NTL)"   = "#E69F00", 
                                  "Predicted Incidence"  = "black")) +
    scale_linetype_manual(values = c("Ecological Risk (Dx)" = "solid", 
                                     "Urban Shield (NTL)"   = "dotted",
                                     "Predicted Incidence"  = "solid")) +
    facet_wrap(~City, ncol = 3) +
    theme_minimal(base_size = 14) +
    labs(title = title_text,
         subtitle = "Smoothed radial profiles showing the decoupling of risk and incidence in cities.",
         x = "Distance from City Center (km)",
         y = "Normalised Intensity",
         color = NULL, linetype = NULL) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  }

# 5. Figure 10: The Urban Paradox (Main Cities) ---------------------------
main_cities <- c("Lagos", "Tamale", "Jos")
plot_data_main <- normalised_df |> 
  filter(City %in% main_cities) |>
  mutate(City = factor(City, levels = main_cities))

p_main <- plot_urban_paradox(plot_data_main, "The Urban Paradox: Hazard vs. Incidence")

ggsave(here(figs_dir, "Fig10_Urban_Paradox_Profiles.png"), p_main, width = 11, height = 5, bg="white")

# 6. Supplementary Figure: All Cities -------------------------------------
p_supp <- plot_urban_paradox(normalised_df |>
                               mutate(City = paste0(City, " (", type, ")")), "Urban Profiles (Focal Cities)")

ggsave(here(figs_dir, "Supp_All_City_Profiles.png"), p_supp, width = 12, height = 12, bg="white")


# --- 7. Expanded quantitative analysis for urban paradox ---------------------
# Load Modern City Locations (Natural Earth)
# Scale 10 gives the highest resolution (includes towns)
wa_cities_sf <- ne_download(scale = 10, type = "populated_places", 
                            category = "cultural", returnclass = "sf") |>
  # Filter spatially to West Africa bounding box (rough crop first for speed)
  st_crop(xmin = -18, xmax = 16, ymin = 4, ymax = 28) |>
  # Filter by Country ISO Codes defined previously
  filter(ADM0_A3 %in% west_africa_countries$iso3c) |>
  select(name = NAME, lat = LATITUDE, lon = LONGITUDE, ADM0_A3)

# Convert to SpatVector for extraction
city_vect <- vect(wa_cities_sf)

# Load Population Raster 
pop_path <- here(proc_dir, "pop_count_2020_05deg.tif")
if(file.exists(pop_path)) {
  Pop_2020 <- terra::rast(pop_path)
} else {
  stop("Population raster required for classification.")
}

# Extract Sum of Population within 10km Radius
city_buffer <- terra::buffer(city_vect, width = 10000) # 10km buffer
city_pop_sums <- terra::extract(Pop_2020, city_buffer, fun = sum, na.rm = TRUE)

all_wa_cities <- wa_cities_sf |>
  mutate(pop_2020 = city_pop_sums[,2]) |>
  st_drop_geometry() |>
  filter(pop_2020 > 20000) |>
  mutate(type_stat = case_when(pop_2020 > 1000000 ~ "Large City",
                               pop_2020 > 300000  ~ "Regional City",
                               pop_2020 > 50000   ~ "Town",
                               TRUE ~ "Small Settlement")) |>
  filter(type_stat != "Small Settlement")

# Check the new distribution
table(all_wa_cities$type_stat)

# Stratified Sample (Up to 5 per Class per Country)
set.seed(123) 
stratified_sample <- all_wa_cities |>
  rename(iso3c = ADM0_A3) |>
  group_by(iso3c, type_stat) |>
  slice_sample(n = 5) |> 
  ungroup() |>
  select(name, lat, lon, pop = pop_2020, type_stat)

# Combine with Focal Cities
final_city_list <- bind_rows(standardised_cities |> 
                               mutate(Dataset = "Focal", 
                                      type_group = case_when(type == "Megacity" ~ "Large City",
                                                             type == "Endemic Town" ~ "Town",
                                                             TRUE ~ type)),
                             stratified_sample |> 
                               mutate(Dataset = "Statistical", type_group = type_stat)) |>
  group_by(name) |>
  slice_head(n = 1) |> 
  ungroup()

all_profiles <- list()
for(i in 1:nrow(final_city_list)) {
  cit <- final_city_list[i,] 
  message(paste("Processing:", cit$name, "(", i, "/", nrow(final_city_list), ")"))
  tryCatch({
    if(!is.na(cit$lon) & !is.na(cit$lat)) {
      all_profiles[[i]] <- get_radial_profile(cit$name, cit$lon, cit$lat, City_Stack)
    }
  }, error = function(e) {
    message(paste("  -> Failed:", cit$name, "-", e$message))
  })
}

# Bind results
expanded_profile_df <- bind_rows(all_profiles) |>
  left_join(final_city_list[, c("name", "type_group", "Dataset")], by = c("City" = "name"))

# Generate Summary Statistics
stats_summary <- expanded_profile_df |>
  filter(Variable == "Incidence") |>
  group_by(City, type_group) |>
  arrange(desc(Mean)) |>
  slice_head(n = 1) |>
  summarise(Peak_Dist = Distance, .groups = "drop")

typology_stats <- stats_summary |>
  group_by(type_group) |>
  summarise(n = n(),
            Median_Peak_Km = median(Peak_Dist),
            IQR_Lower = quantile(Peak_Dist, 0.25),
            IQR_Upper = quantile(Peak_Dist, 0.75)) |>
  mutate(Text = paste0("n = ", n, ", median = ", round(Median_Peak_Km, 0), 
                  " km, IQR = ", round(IQR_Lower, 0), "–", round(IQR_Upper, 0), " km")) |>
  arrange(desc(Median_Peak_Km))

kw_test <- kruskal.test(Peak_Dist ~ type_group, data = stats_summary)
