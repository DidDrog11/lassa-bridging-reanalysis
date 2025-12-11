################################################################################
## SCRIPT: 1_descriptive_analysis.R
## PURPOSE: Generate Table 1 (Data Summary) and Figure 1 (Data Distribution Map).
## INPUTS:  Cleaned Rodent Data, Pathogen Training Data, Human Calibration Data,
##          Validation Case Counts, NTL Raster (for Urban/Rural calc).
## OUTPUTS: results/figures/Fig1_Data_Distribution.png
##          results/tables/Table1_Data_Summary.csv
################################################################################

source(here::here("packages.R"))

# Directories
maps_dir <- here("results", "maps")
figs_dir <- here("results", "figures")
tabs_dir <- here("results", "tables")
proc_dir <- here("data", "processed")
raw_dir  <- here("data", "raw")

if(!dir.exists(figs_dir)) dir.create(figs_dir)
if(!dir.exists(tabs_dir)) dir.create(tabs_dir)

# 1. Load Data ------------------------------------------------------------
# A. Geospatial Context (West Africa)
west_africa_countries <- countrycode::codelist |>
  filter(str_detect(region23, "Western Africa")) |>
  select(country = country.name.en, iso3c, iso2c) |>
  filter(!str_detect(iso3c, "CPV|SHN"))
west_africa_countries$country[west_africa_countries$country == "Côte d’Ivoire"] <- "Ivory Coast"

wa_vect <- ne_countries(scale = 10, country = west_africa_countries$country, returnclass = "sv")

country_labels <- terra::centroids(wa_vect, inside=TRUE) |>
  as.data.frame(geom = "XY") |>
  rename(lon = x, lat = y, label = admin)

africa_sf <- ne_countries(scale = 110, continent = "Africa", returnclass = "sf") |>
  mutate(Region = ifelse(admin %in% west_africa_countries$country, "West Africa", "Rest of Africa")) |>
  select(adm0_a3, Region)

# B. Host Data (Rodent Occurrences)
# From clean_data.R output
rodent_df <- read_rds(here(proc_dir, "rodent_data_mapped.rds")) |>
  mutate(Plot_Group = case_when(
    species == "Mastomys natalensis" ~ "Primary Reservoir (M. natalensis)",
    species %in% c("Rattus rattus", "Mus musculus") ~ "Invasive (Rattus/Mus)",
    TRUE ~ "Native Community"
  )) |>
  # Set Factor Levels for Plotting Order (Native -> Invasive -> Reservoir)
  mutate(Plot_Group = factor(Plot_Group, levels = c("Native Community", 
                                                    "Invasive (Rattus/Mus)", 
                                                    "Primary Reservoir (M. natalensis)"))) |>
  mutate(Source_Category = case_when(
    source %in% c("ArHa_PA", "WA_Rodents_PA") ~ "Systematic Surveys (ArHa/WA)",
    source %in% c("UrbanLit_HQ", "UrbanLit_PO") ~ "Targeted Urban Search",
    source %in% c("GBIF_PO", "Original_PO", "Opportunistic_PO") ~ "Opportunistic / GBIF",
    TRUE ~ "Other"
  ))

iucn_path <- here(raw_dir, "iucn", "data_0.shp")
if(file.exists(iucn_path)) {
  iucn_vect <- vect(iucn_path)
} else {
  warning("IUCN Shapefile not found.")
  iucn_vect <- NULL
}

# C. Pathogen Training Data (Rodent Testing)
pathogen_rodent <- readRDS(here(proc_dir, "pathogen_training_data.rds")) |>
  mutate(assay_source = case_when(source == "SL_pathogen" ~ "Serology",
                                  source == "Basinski_pathogen" ~ "PCR",
                                  TRUE ~ NA))
# Pathogen raw data for assay origin
arha_path <- here(raw_dir, "Project_ArHa_database.rds")
if (!file.exists(arha_path)) stop("Project ArHa RDS not found.")
arha_db <- readRDS(arha_path)
lasv_rodent_arha <- arha_db$pathogen |>
  left_join(arha_db$host, by = "host_record_id") |>
  filter(host_species == "Mastomys natalensis",
         pathogen_species_cleaned == "Mammarenavirus lassaense",
         !is.na(latitude), !is.na(longitude),
         number_tested > 0) |>
  # Aggregate to Site Level
  group_by(latitude, longitude, assay) |>
  summarise(n_tested = sum(number_tested, na.rm = TRUE),
            n_pos = sum(number_positive, na.rm = TRUE),
            .groups = "drop") |>
  mutate(outcome = case_when(n_pos > 0 ~ 1, # Define Outcome: 1 = Found, 0 = Not Found (if effort sufficient)
                             n_tested >= 1 ~ 0,
                             TRUE ~ NA_real_),
         source = "ArHa_pathogen") |> # Ambiguous
  filter(!is.na(outcome)) |>
  select(latitude, longitude, assay, source)

pathogen_rodent <- left_join(pathogen_rodent, lasv_rodent_arha) |>
  mutate(assay = coalesce(assay, assay_source)) |>
  select(-assay_source)

# D. Human Calibration Data (Serosurveys + Anchors)
human_df <- read.csv(here(proc_dir, "human_calibration_data.csv"))

# E. Validation Data (All 4 Sets)
# 1. Binary LGA (NCDC) - High Res Nigeria
val_binary_lga <- read.csv(here(proc_dir, "validation_binary_lga.csv")) |>
  group_by(LGA_Code) |>
  summarise(Status = max(Status, na.rm = TRUE))

# 2. Case Counts Admin 2 (Moore et al) - Regional
val_cases_adm2 <- read.csv(here(proc_dir, "validation_cases_admin2.csv"))

# 3. Case Counts Admin 1 (General)
val_cases_adm1 <- read.csv(here(proc_dir, "validation_cases_admin1.csv"))

# 4. Nigeria Admin 1 (NCDC Aggregate)
val_cases_nga1 <- read.csv(here(proc_dir, "validation_cases_nigeria_adm1.csv"))

# F. NTL Raster
# Assuming final stack exists. If not, try loading raw or skipping.
stack_path <- here(proc_dir, "final_predictor_stack_scaled.tif")
if(file.exists(stack_path)) {
  NTL_rast <- terra::rast(stack_path)[["NTL"]]
} else {
  warning("NTL Raster not found. Urban/Rural counts will be NA.")
  NTL_rast <- NULL
}

# G. Pop Class Raster
degurba_path <- here(proc_dir, "pop_class_dijkstra_05deg.tif") 

if(file.exists(degurba_path)) {
  DegUrba_rast <- terra::rast(degurba_path)
  message("Loaded Dijkstra Urban Classification Raster.")
} else {
  warning("Dijkstra Raster not found. Urban counts will be NA.")
  DegUrba_rast <- NULL
}

# Target cities
standardised_cities <- read_rds(here(proc_dir, "standardised_cities.rds"))

# 2. Table 1: Data Summary ------------------------------------------------
classify_context <- function(lat, lon) {
  if(is.null(DegUrba_rast)) return(c(Rural=NA, Town=NA, City=NA))
  
  pts <- terra::vect(cbind(lon, lat), crs="EPSG:4326")
  vals <- terra::extract(DegUrba_rast, pts, ID=FALSE)[,1]
  
  # Your definition: 1=Rural, 2=Town, 3=City
  n_rur  <- sum(vals == 1, na.rm=TRUE)
  n_town <- sum(vals == 2, na.rm=TRUE)
  n_city <- sum(vals == 3, na.rm=TRUE)
  
  return(c(Rural=n_rur, Town=n_town, City=n_city))
}

# A. Summarise Hosts
if(!is.null(DegUrba_rast)) {
  pts_all <- terra::vect(cbind(rodent_df$lon, rodent_df$lat), crs="EPSG:4326")
  vals_all <- terra::extract(DegUrba_rast, pts_all, ID=FALSE)[,1]
  rodent_df$Context <- case_when(
    vals_all == 3 ~ "City",
    vals_all == 2 ~ "Town",
    vals_all == 1 ~ "Rural",
    TRUE ~ "Unknown"
  )
} else {
  rodent_df$Context <- "Unknown"
}

site_context_summary <- rodent_df |>
  distinct(Source_Category, lat, lon, Context) |> # Collapse to unique sites
  group_by(Source_Category) |>
  summarise(
    Rural_Sites = sum(Context == "Rural", na.rm=TRUE),
    Town_Sites  = sum(Context == "Town", na.rm=TRUE),
    City_Sites  = sum(Context == "City", na.rm=TRUE),
    .groups = "drop"
  )

# Generate Stratified Table
tab_host <- rodent_df |>
  group_by(Source_Category) |>
  summarise(
    Target = "Small Mammal Community",
    N_Sites = n_distinct(paste(lat, lon)),
    # N_Samples = Total Species-Site Detections (Occupancy Data Points)
    N_Samples = n(), 
    `N_Mnat (1 sp.)` = sum(species == "Mastomys natalensis"),
    `N_Invasive (2 spp.)` = sum(species %in% c("Rattus rattus", "Mus musculus")),
    `N_Native (5 spp.)` = sum(!species %in% c("Mastomys natalensis", "Rattus rattus", "Mus musculus")),
    Period = paste0(min(year, na.rm = TRUE), "-", max(year, na.rm = TRUE))
  ) |>
  left_join(site_context_summary, by = "Source_Category") |>
  arrange(desc(City_Sites))

# B. Summarise Pathogen (Rodent)
if(!is.null(DegUrba_rast)) {
  pts_path <- terra::vect(cbind(pathogen_rodent$longitude, pathogen_rodent$latitude), crs = "EPSG:4326")
  vals_path <- terra::extract(DegUrba_rast, pts_path, ID = FALSE)[,1]
  pathogen_rodent$Context <- case_when(
    vals_path == 3 ~ "City",
    vals_path == 2 ~ "Town",
    vals_path == 1 ~ "Rural",
    TRUE ~ "Unknown"
  )
} else {
  pathogen_rodent$Context <- "Unknown"
}

tab_path <- pathogen_rodent |>
  group_by(assay) |>
  summarise(
    Target = "M. natalensis",
    N_Sites = n_distinct(paste(latitude, longitude)),
    N_Samples = sum(n_tested),
    N_Positive = sum(n_pos),
    Prevalence = paste0(round(sum(n_pos)/sum(n_tested)*100, 1), "%"),
    Rural_Sites = sum(Context == "Rural", na.rm = TRUE),
    Town_Sites  = sum(Context == "Town", na.rm = TRUE),
    City_Sites  = sum(Context == "City", na.rm = TRUE),
    `N_Mnat (1 sp.)` =  sum(n_tested),
    `N_Invasive (2 spp.)` = NA,
    `N_Native (5 spp.)` = NA,
    Period = "1972-2022" ) |>
  rename(Source_Category = assay)

# C. Summarise Human
real_human <- human_df |>
  filter(str_detect(source, "pathogen"))
syn_human <- human_df |>
  filter(str_detect(source, "Synthetic"))

# Calculate Context
ctx_hum <- classify_context(real_human$lat, real_human$lon)
ctx_syn <- classify_context(syn_human$lat, syn_human$lon)

total_rural <- ctx_hum["Rural"] + ctx_syn["Rural"]
total_town  <- ctx_hum["Town"]  + ctx_syn["Town"]
total_city  <- ctx_hum["City"]  + ctx_syn["City"]

tab_human <- human_df |>
  summarise(
    Source_Category = "Human Seroprevalence",
    Target = "Community (IgG)",
    N_Sites = n(), 
    N_Samples = sum(n_test),
    N_Positive = sum(n_pos),
    Prevalence = paste0(round(sum(n_pos)/sum(n_test)*100, 1), "%"),
    Rural_Sites = total_rural,
    Town_Sites  = total_town,
    City_Sites  = total_city, 
    N_Mnat = NA, N_Invasive = NA, N_Native = NA,
    Period = paste0(min(real_human$year, na.rm = TRUE), "-", max(real_human$year, na.rm = TRUE))
  )

# Validation data summary

tab_val_nga <- data.frame(
  Source_Category = "NCDC Surveillance (Nigeria)",
  Target = "Confirmed Cases (Binary)",
  Period = "2018-2025",
  N_Sites = n_distinct(val_binary_lga$LGA_Code),
  N_Samples = NA, N_Positive = sum(val_binary_lga$Status), Prevalence = NA,
  City_Sites = NA, Town_Sites = NA, Rural_Sites = NA,
  N_Mnat = NA, N_Invasive = NA, N_Native = NA
)

tab_val_moore <- data.frame(
  Source_Category = "Regional Surveillance (Moore et al.)",
  Target = "Confirmed Cases (Count)",
  Period = "2012-2022",
  N_Sites = n_distinct(val_cases_adm2$Location_Code),
  N_Samples = NA, N_Positive = sum(val_cases_adm2$Annual_Cases > 0), Prevalence = NA,
  City_Sites = NA, Town_Sites = NA, Rural_Sites = NA,
  N_Mnat = NA, N_Invasive = NA, N_Native = NA
)

colnames(tab_human)[which(names(tab_human) == "N_Mnat")] <- "N_Mnat (1 sp.)"
colnames(tab_human)[which(names(tab_human) == "N_Invasive")] <- "N_Invasive (2 spp.)"
colnames(tab_human)[which(names(tab_human) == "N_Native")] <- "N_Native (5 spp.)"

colnames(tab_val_nga)[which(names(tab_val_nga) == "N_Mnat")] <- "N_Mnat (1 sp.)"
colnames(tab_val_nga)[which(names(tab_val_nga) == "N_Invasive")] <- "N_Invasive (2 spp.)"
colnames(tab_val_nga)[which(names(tab_val_nga) == "N_Native")] <- "N_Native (5 spp.)"

colnames(tab_val_moore)[which(names(tab_val_moore) == "N_Mnat")] <- "N_Mnat (1 sp.)"
colnames(tab_val_moore)[which(names(tab_val_moore) == "N_Invasive")] <- "N_Invasive (2 spp.)"
colnames(tab_val_moore)[which(names(tab_val_moore) == "N_Native")] <- "N_Native (5 spp.)"

# Combine All
table_1 <- bind_rows(tab_host, tab_path, tab_human, tab_val_nga, tab_val_moore) |>
  select(Source_Category, Target, Period, N_Sites, N_Positive, 
         City_Sites, Town_Sites, Rural_Sites, 
         `N_Mnat (1 sp.)`, `N_Invasive (2 spp.)`, `N_Native (5 spp.)`)

print(table_1)
write.csv(table_1, here(tabs_dir, "Table1_Data_Summary.csv"), row.names = FALSE)

# 3. Figure 1: The Three-Panel Map ----------------------------------------
# --- Panel A: Host Sampling ---
rodent_plot_df <- rodent_df |>
  filter(source != "Opportunistic_PO") |>
  arrange(Plot_Group)

p_host <- ggplot() +
  geom_spatvector(data = wa_vect, fill = "gray97", colour = "gray60", linewidth = 0.8) +
  {if(!is.null(iucn_vect)) geom_spatvector(data = crop(iucn_vect, wa_vect), fill = NA, colour = "black", linewidth = 0.6, linetype = "dashed")} +
  geom_point(data = rodent_plot_df, 
             aes(x = lon, y = lat, 
                 colour = Plot_Group, size = Plot_Group, alpha = Plot_Group),
             position = position_jitter(width = 0.05, height = 0.05)) +
  scale_colour_manual(values = c("Primary Reservoir (M. natalensis)" = "#D55E00", 
                                 "Invasive (Rattus/Mus)" = "#0072B2",              
                                 "Native Community" = "gray50")) +
  scale_size_manual(values = c("Primary Reservoir (M. natalensis)" = 1.2, 
                               "Invasive (Rattus/Mus)" = 1.2,              
                               "Native Community" = 1)) +
  scale_alpha_manual(values = c("Primary Reservoir (M. natalensis)" = 0.7, 
                                "Invasive (Rattus/Mus)" = 0.6,              
                                "Native Community" = 0.4)) +
  theme_minimal() +
  labs(tag = "A", subtitle = "Host Community Sampling (8 Species)", colour = NULL, size = NULL, alpha = NULL, x = NULL, y = NULL) +
  theme(legend.position = c(0.65, 0.85), 
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        axis.text = element_blank())

# --- Panel B: Pathogen Sampling ---
path_pts <- bind_rows(pathogen_rodent |>
                        select(lon = longitude, lat = latitude) |>
                        mutate(Type = "Rodent Testing"),
                      human_df |>
                        filter(source != "Synthetic_Constraint") |>
                        select(lon, lat) |>
                        mutate(Type = "Human Serosurvey"))

# Add Anchors specifically to highlight them
anchor_pts <- human_df |>
  filter(source == "Synthetic_Constraint") |>
  select(lon, lat) |>
  mutate(Type = "Urban Anchor")

p_pathogen <- ggplot() +
  geom_spatvector(data = wa_vect, fill = "gray90", color = "black", linewidth = 0.8) +
  geom_text_repel(data = country_labels, aes(x = lon, y = lat, label = label),
                  size = 3, fontface = "bold", colour = "gray20",
                  bg.color = "white", bg.r = 0.15,
                  seed = 123, box.padding = 0.5) +
  geom_point(data = path_pts, aes(x = lon, y = lat, colour = Type, shape = Type), 
             size = 1, stroke = 0.8, alpha = 0.8,
             position = position_jitter(width = 0.05, height = 0.05)) +
  geom_point(data = anchor_pts, aes(x = lon, y = lat, shape = Type), 
             colour = "black", size = 3, stroke = 1) +
  scale_colour_manual(values = c("Rodent Testing" = "firebrick", 
                                 "Human Serosurvey" = "gold3", 
                                 "Urban Anchor" = "black")) +
  scale_shape_manual(values = c("Rodent Testing" = 16, 
                                "Human Serosurvey" = 17, 
                                "Urban Anchor" = 8)) +
  theme_minimal() +
  labs(tag = "B", subtitle = "LASV Prevalence Data", colour = NULL, shape = NULL, x = NULL, y = NULL) +
  theme(legend.position = c(0.65, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        axis.text = element_blank())

# --- Panel C: Validation (Cases) ---
adm_2_files <- list.files(path = here(raw_dir, "gadm"), pattern = "_2_", full.names = TRUE)

if(length(adm_2_files) > 0) {
  adm_2_vect <- do.call(rbind, lapply(adm_2_files, vect)) 
  adm_2_cases <- adm_2_vect |> 
    left_join(val_cases_adm2, by = c("GID_2" = "Location_Code")) |>
    filter(!is.na(Annual_Cases) & Annual_Cases > 0 & !str_detect(GID_0, "NGA"))
  nga_adm2 <- adm_2_vect[adm_2_vect$GID_0 == "NGA", ]
  nga_binary <- nga_adm2 |>
    left_join(val_binary_lga, by = c("GID_2" = "LGA_Code")) |>
    filter(Status == 1) # Only show presence
}

relevant_isos <- unique(c(adm_2_cases$GID_0, nga_binary$GID_0))

adm_2_context <- adm_2_vect |>
  filter(GID_0 %in% relevant_isos)

standardised_cities_sf <- vect(standardised_cities, geom = c("lon", "lat"), crs = "EPSG:4326") |>
  st_as_sf()

p_cases <- ggplot() +
  geom_spatvector(data = wa_vect, fill = "gray97", colour = "black") +
  geom_spatvector(data = adm_2_context, fill = "gray70", colour = "black", linewidth = 0.05) +
  geom_spatvector(data = wa_vect, fill = NA, colour = "black", linewidth = 0.8) +
  {if(exists("nga_binary")) geom_spatvector(data = nga_binary, fill = "firebrick", colour = "black")} +
  {if(exists("adm_2_cases")) geom_spatvector(data = adm_2_cases, aes(fill = log1p(Annual_Cases)), colour = NA)} +
  geom_label_repel(data = standardised_cities_sf, aes(label = name, geometry = geometry),
                   stat = "sf_coordinates",
                   min.segment.length = 0) +
  scale_fill_viridis_c(option = "magma", name = "Log Cases") +
  theme_minimal() +
  labs(tag = "C", subtitle = "Validation: NCDC Confirmed (Red) & Regional Reports", x = NULL, y = NULL) +
  theme(legend.position = c(0.65, 0.85), legend.background = element_rect(fill = alpha("white", 0.8), colour = NA), axis.text = element_blank())

p_inset <- ggplot(africa_sf) +
  geom_sf(aes(fill = Region), colour = "white", size = 0.1) +
  scale_fill_manual(values = c("Rest of Africa" = "gray85", "West Africa" = "#D55E00")) +
  theme_void() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

# --- Assemble ---
# Combine A and B side-by-side
top_row <- p_host | p_pathogen

# Place inset into C (top-left empty space)
bottom_row <- p_cases +
  inset_element(p_inset, left = 0.02, bottom = 0.6, right = 0.3, top = 0.98, align_to = 'plot')

# Final layout
final_plot <- top_row / bottom_row +
  plot_layout(heights = c(1, 1.3)) +
  theme(plot.tag = element_text(face = "bold", size = 14))

# Save (Wider width for side-by-side layout)
ggsave(here(figs_dir, "Fig1_Data_Distribution.png"), final_plot, width = 12, height = 12, dpi = 300, bg = "white")
