################################################################################
## SCRIPT: 4_risk_and_burden_interpretation.R
## PURPOSE: 1. Visualize the "Socio-economic Shield" (Calibration Curves).
##          2. Map Final Incidence (Cases per pixel).
##          3. Generate Country-level Burden Table.
##          4. Validate against NCDC Data (Silent Districts).
## INPUTS:  Dm/Dl rasters, Incidence_Full_Custom_Lambda0.03.tif,
##          human_calibration_data.csv, method_comparison_summary.rds
## OUTPUTS: Fig7_Calibration_Shield.png, Fig8_Incidence_Map.png,
##          Fig9_Silent_Districts.png, Table2_Country_Burden.csv
################################################################################

source(here::here("packages.R"))

# Directories
maps_dir <- here("results", "maps")
figs_dir <- here("results", "figures")
tabs_dir <- here("results", "tables")
proc_dir <- here("data", "processed")
data_dir <- here("data", "raw")

if(!dir.exists(figs_dir)) dir.create(figs_dir)

# A. Rasters (Custom & Basinski)
# Host
Dm_custom <- terra::rast(here(maps_dir, "Dm_Mastomys_natalensis.tif"))
bas_dm_path <- here(data_dir, "basinski_comparison", "host_raster.tif")
Dm_basinski <- terra::rast(bas_dm_path)

# Pathogen
Dl_custom <- terra::rast(here(maps_dir, "Dl_Pathogen_Layer_Ensemble.tif"))
bas_dl_path <- here(data_dir, "basinski_comparison", "pathogen_raster.tif")
Dl_basinski <- terra::rast(bas_dl_path)

# B. IUCN Mask (The Biological Filter)
iucn_path <- here(data_dir, "iucn", "data_0.shp")

if(file.exists(iucn_path)) {
  iucn_vect <- vect(iucn_path)
  iucn_mask <- terra::rasterize(iucn_vect, Dm_custom, field=1, background=NA)
  
  if(!terra::compareGeom(Dm_custom, Dm_basinski, stopOnError=FALSE)) {
    Dm_basinski <- terra::resample(Dm_basinski, Dm_custom, method="bilinear")
    Dl_basinski <- terra::resample(Dl_basinski, Dm_custom, method="bilinear")
  }
  
  Dm_custom   <- terra::mask(Dm_custom, iucn_mask)
  Dl_custom   <- terra::mask(Dl_custom, iucn_mask)
  Dm_basinski <- terra::mask(Dm_basinski, iucn_mask)
  Dl_basinski <- terra::mask(Dl_basinski, iucn_mask)
  
  message("IUCN Mask applied to all rasters.")
} else {
  warning("IUCN Shapefile not found. Maps will be unmasked!")
}

# Calculate Dx
Dx_custom   <- Dm_custom * Dl_custom
Dx_basinski <- Dm_basinski * Dl_basinski
names(Dx_custom) <- "Dx"
names(Dx_basinski) <- "Dx"

# Load Final Incidence Raster
inc_path <- here(maps_dir, "Incidence_Full_Custom_Lambda0.03.tif")
if(file.exists(inc_path)) {
  Incidence_Custom <- terra::rast(inc_path)
  # Ensure incidence is also masked
  if(exists("iucn_mask")) Incidence_Custom <- terra::mask(Incidence_Custom, iucn_mask)
} else {
  stop("Incidence Raster not found. Please run 'calculate_risk_and_validate.R' first.")
}

# C. Context Layers
west_africa_countries <- countrycode::codelist |>
  filter(str_detect(region23, "Western Africa")) |>
  select(country = country.name.en, iso3c) |>
  filter(!str_detect(iso3c, "CPV|SHN"))
west_africa_countries$country[west_africa_countries$country == "Côte d’Ivoire"] <- "Ivory Coast"
wa_vect <- ne_countries(scale = 10, country = west_africa_countries$country, returnclass = "sv")

# D. Human Data
human_df <- read.csv(here(proc_dir, "human_calibration_data.csv")) |>
  filter(n_test > 0)
human_vect <- vect(human_df, geom = c("lon", "lat"), crs = "EPSG:4326")

stack_path <- here(proc_dir, "final_predictor_stack_scaled.tif")
if(file.exists(stack_path)) {
  NTL_rast <- terra::rast(stack_path)[["NTL"]]
  human_df$Dx_Val  <- terra::extract(Dx_custom, human_vect, ID=FALSE)[,1]
  human_df$NTL_Val <- terra::extract(NTL_rast, human_vect, ID=FALSE)[,1]
  human_df <- na.omit(human_df |>
                        select(-id))
}
# update after removing NAs
human_vect <- vect(human_df, geom = c("lon", "lat"), crs = "EPSG:4326")

# E. Validation Data
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

# 2. Figure 7: The Socio-economic Shield (Calibration) --------------------
degurba_path <- here(proc_dir, "pop_class_dijkstra_05deg.tif")
DegUrba_rast <- terra::rast(degurba_path)

# Extract Urban Class for points
human_df$Urban_Class <- terra::extract(DegUrba_rast, human_vect, ID=FALSE)[,1]

# Define Representative NTL Values based on Dijkstra Classes
ntl_summary <- human_df |>
  group_by(Urban_Class) |>
  summarise(Median_NTL = median(NTL_Val, na.rm=TRUE))
rural_ntl <- ntl_summary$Median_NTL[ntl_summary$Urban_Class == 1]
urban_ntl <- ntl_summary$Median_NTL[ntl_summary$Urban_Class == 3]

# Fit the GLM 1000 times to get robust Confidence Intervals
n_boot <- 1000
boot_preds <- list()

set.seed(123)

# Create prediction grid
pred_grid <- expand.grid(
  Dx_Val = seq(0, 1, length.out = 100),
  NTL_Val = c(rural_ntl, urban_ntl)
)
pred_grid$Context <- ifelse(pred_grid$NTL_Val == rural_ntl, "Rural", "Urban")

for(i in 1:n_boot) {
  # Resample data with replacement
  df_boot <- human_df[sample(nrow(human_df), replace = TRUE), ]
  
  # Fit GLM (Quasibinomial to handle the dispersion mechanics, though resampling handles the variance)
  # We use a simple additive model. Interaction (Dx * NTL) could be tested but might be unstable.
  m_boot <- tryCatch(
    glm(cbind(n_pos, n_test - n_pos) ~ Dx_Val + NTL_Val, 
        family = binomial(link = "logit"), # Binomial fine for coefficients in bootstrap
        data = df_boot),
    error = function(e) NULL
  )
  
  if(!is.null(m_boot)) {
    # Predict
    p <- predict(m_boot, newdata = pred_grid, type = "response")
    boot_preds[[i]] <- data.frame(
      Iter = i,
      Context = pred_grid$Context,
      Dx_Val = pred_grid$Dx_Val,
      Predicted = p
    )
  }
}

boot_df <- bind_rows(boot_preds) |>
  group_by(Context, Dx_Val) |>
  summarise(
    Mean_Pred = mean(Predicted),
    Lower_CI = quantile(Predicted, 0.025),
    Upper_CI = quantile(Predicted, 0.975),
    .groups = "drop"
  )

binned_data <- human_df |>
  mutate(Context = case_when(Urban_Class == 1 ~ "Rural", Urban_Class == 3 ~ "Urban", TRUE ~ NA_character_)) |>
  filter(!is.na(Context)) |>
  mutate(Dx_Bin = cut(Dx_Val, breaks = 10)) |>
  group_by(Context, Dx_Bin) |>
  summarise(
    Mean_Dx = mean(Dx_Val),
    Observed_Prev = sum(n_pos) / sum(n_test),
    Total_N = sum(n_test),
    .groups = "drop"
  )

p_shield <- ggplot() +
  geom_ribbon(data = boot_df, aes(x = Dx_Val, ymin = Lower_CI, ymax = Upper_CI, fill = Context), alpha = 0.2) +
  geom_line(data = boot_df, aes(x = Dx_Val, y = Mean_Pred, colour = Context, linetype = Context), size = 1.2) +
  geom_point(data = binned_data, aes(x = Mean_Dx, y = Observed_Prev, color = Context, size = Total_N), alpha = 0.9) +
  scale_colour_manual(values = c("Rural" = "#D55E00", "Urban" = "#0072B2")) +
  scale_fill_manual(values = c("Rural" = "#D55E00", "Urban" = "#0072B2")) +
  scale_linetype_manual(values = c("Rural" = "solid", "Urban" = "dashed")) +
  scale_size_area(max_size = 6) + 
  coord_cartesian(ylim = c(0, 0.6)) + 
  theme_minimal(base_size = 14) +
  labs(title = "The Socio-economic Shield",
       x = "Ecological Risk Index (Dx)", 
       y = "Predicted Seroprevalence",
       size = "Sample Size", color = NULL, fill = NULL, linetype = NULL) +
  theme(legend.position = "bottom")

ggsave(here(figs_dir, "Fig7_Calibration_Shield.png"), p_shield, width = 8, height = 6, bg="white")

# 3. Figure 8: The Incidence Map (Main Result) ----------------------------
p_inc <- ggplot() +
  geom_spatraster(data = Incidence_Custom) +
  geom_spatvector(data = wa_vect, fill = NA, color = "black", linewidth = 0.2) +
  scale_fill_viridis_c(option = "turbo", name = "Annual Infections\n(per 0.05° Cell)", 
                       trans = "log1p", breaks = c(0, 10, 100, 1000, 3000), 
                       labels = c("0", "10", "100", "1,000", "3,000"),
                       na.value = "transparent") +
  theme_minimal() +
  labs(title = "Predicted Annual LASV Infections",
       x = NULL, y = NULL) +
  theme(legend.position = "right", legend.key.height = unit(1.5, "cm"))

ggsave(here(figs_dir, "Fig8_Incidence_Map.png"), p_inc, width = 8, height = 6, bg="white")

# 4. Table 2: Burden by Country -------------------------------------------
adm0_rast <- terra::rasterize(wa_vect |>
                                select(sov_a3, admin), Incidence_Custom, field = "admin")

pop_path <- here(proc_dir, "pop_count_2020_05deg.tif")
if(file.exists(pop_path)) {
  Pop_rast <- terra::rast(pop_path)
} else { stop("Population raster missing.") }

# Calculate Zonal Sums
cases_by_country <- terra::zonal(Incidence_Custom, adm0_rast, fun = "sum", na.rm = TRUE) |>
  rename(Country = admin, Total_Cases = Omega_Full_Custom) |>
  drop_na()

pop_by_country <- terra::zonal(Pop_rast, adm0_rast, fun = "sum", na.rm = TRUE) |>
  rename(Country = admin, Total_Pop = Pop_Count_Raw) |>
  drop_na()

# Format Table
table_2 <- cases_by_country |>
  left_join(pop_by_country, by = "Country") |>
  mutate(`1000s_of_Infections` = round(Total_Cases / 1000, 1),
         Rate = round((Total_Cases / Total_Pop) * 1000, 1)) |>
  arrange(desc(Total_Cases)) |>
  select(Country, `1000s_of_Infections`, Rate)

print(table_2)
write.csv(table_2, here(tabs_dir, "Table2_Country_Burden.csv"), row.names = FALSE)


# 5. Figure 9: Surveillance Gaps (Split Panels) ---------------------------
# --- PANEL A: NIGERIA (NCDC Validation) ---
nga_adm2_path <- list.files(here(data_dir, "gadm"), pattern = "NGA_2", full.names = TRUE)[1]

if(!is.na(nga_adm2_path)) {
  nga_vect <- vect(nga_adm2_path)
  
  # Extract Incidence
  nga_vect$Mean_Inc <- terra::extract(Incidence_Custom, nga_vect, fun = mean, na.rm = TRUE)[,2]
  
  # Join NCDC Data
  nga_df <- values(nga_vect) |>
    left_join(val_binary_lga, by = c("GID_2" = "LGA_Code")) |>
    mutate(Status = tidyr::replace_na(Status, 0)) |>
    mutate(
      Inc_Quantile = ntile(Mean_Inc, 4),
      Category = case_when(
        Status == 1 ~ "Confirmed Endemic",
        Status == 0 & Inc_Quantile == 4 ~ "Silent (High Risk)",
        Status == 0 & Inc_Quantile %in% c(2,3) ~ "Silent (Medium Risk)",
        Status == 0 & Inc_Quantile == 1 ~ "Silent (Low Risk)",
        TRUE ~ "Other"
      )
    ) |>
    mutate(Category = factor(Category, levels = c("Confirmed Endemic", "Silent (High Risk)", "Silent (Medium Risk)", "Silent (Low Risk)")))
  
  values(nga_vect) <- nga_df
  
  # Plot Nigeria
  p_nga <- ggplot() +
    geom_spatvector(data = nga_vect, aes(fill = Category), color = "white", linewidth = 0.05) +
    scale_fill_manual(values = c("Confirmed Endemic" = "firebrick",
                                 "Silent (High Risk)" = "#E69F00", 
                                 "Silent (Medium Risk)" = "#F0E442",
                                 "Silent (Low Risk)" = "gray90")) +
    theme_void() +
    labs(title = "A. Nigeria", subtitle = "Risk stratification of LGAs with zero reported cases.") +
    theme(legend.position = "none") # Collect legend later
  
  # Stats Check (Nigeria)
  auc_nga <- pROC::auc(nga_df$Status, nga_df$Mean_Inc)
  message(paste("Nigeria Prediction AUC:", round(auc_nga, 3)))
}

# --- PANEL B: MANO RIVER REGION (Moore et al. Validation) ---
# Load Guinea, Sierra Leone, Liberia
mru_iso <- c("GIN", "SLE", "LBR")
adm_2_files <- list.files(path = here(data_dir, "gadm"), pattern = "_2_", full.names = TRUE)
# Filter for MRU files
mru_files <- adm_2_files[grepl(paste(mru_iso, collapse="|"), adm_2_files)]

if(length(mru_files) > 0) {
  mru_vect <- do.call(rbind, lapply(mru_files, vect))
  
  # Extract Incidence
  mru_vect$Mean_Inc <- terra::extract(Incidence_Custom, mru_vect, fun = mean, na.rm = TRUE)[,2]
  
  # Join Moore Data
  mru_df <- values(mru_vect) |>
    left_join(val_cases_adm2, by = c("GID_2" = "Location_Code")) |>
    mutate(Annual_Cases = tidyr::replace_na(Annual_Cases, 0)) |>
    mutate(
      Status = ifelse(Annual_Cases > 0, 1, 0),
      Inc_Quantile = ntile(Mean_Inc, 4),
      Category = case_when(
        Status == 1 ~ "Confirmed Endemic",
        Status == 0 & Inc_Quantile == 4 ~ "Silent (High Risk)",
        Status == 0 & Inc_Quantile %in% c(2,3) ~ "Silent (Medium Risk)",
        Status == 0 & Inc_Quantile == 1 ~ "Silent (Low Risk)",
        TRUE ~ "Other"
      )
    ) |>
    mutate(Category = factor(Category, levels = c("Confirmed Endemic", "Silent (High Risk)", "Silent (Medium Risk)", "Silent (Low Risk)")))
  
  mru_vect <- left_join(mru_vect, mru_df, by = c("GID_2"))
  
  # Plot MRU
  p_mru <- ggplot() +
    geom_spatvector(data = mru_vect, aes(fill = Category), colour = "white", linewidth = 0.05) +
    scale_fill_manual(values = c("Confirmed Endemic" = "firebrick",
                                 "Silent (High Risk)" = "#E69F00", 
                                 "Silent (Medium Risk)" = "#F0E442",
                                 "Silent (Low Risk)" = "gray90")) +
    geom_spatvector(data = wa_vect |>
                      filter(sov_a3 %in% mru_vect$GID_0.x), fill = "transparent", colour = "black", linewidth = 0.01) +
    theme_void() +
    labs(title = "B. Mano River Region", subtitle = "Validation against regional case counts.") +
    theme(legend.position = "right")
  
  # Stats Check (MRU)
  auc_mru <- pROC::auc(mru_df$Status, mru_df$Mean_Inc)
  message(paste("Mano River Prediction AUC:", round(auc_mru, 3)))
}

# --- COMBINE AND SAVE ---
if(exists("p_nga") & exists("p_mru")) {
  # Layout: Nigeria (A) | MRU (B)
  # We use 'collect' to share the legend from B
  fig9 <- (p_nga | p_mru) + 
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Surveillance Gaps in West Africa",
      subtitle = "Identifying districts with high predicted incidence but no reported cases."
    ) & theme(legend.position = "bottom")
  
  ggsave(here(figs_dir, "Fig9_Silent_Districts.png"), fig9, width = 12, height = 8, bg="white")
}

# 6. Load incidence burden ---------------------------------------------------
summary_path <- here(tabs_dir, "method_comparison_summary.rds")

if(file.exists(summary_path)) {
  burden_summary <- readRDS(summary_path)
  
  # Extract the specific values for the "Full_Custom" (Shielded) model
  stats <- burden_summary |>
    dplyr::filter(Scenario == "Full_Custom", Model_Type == "Shielded") |>
    dplyr::select(Lambda, Total_Infections)
  
  # Format numbers for text (Millions with 1 decimal)
  headline <- round(stats$Total_Infections[stats$Lambda == 0.03] / 1e6, 1)
  lower    <- round(stats$Total_Infections[stats$Lambda == 0] / 1e6, 1)
  upper    <- round(stats$Total_Infections[stats$Lambda == 0.064] / 1e6, 1)
  
  message("------------------------------------------------------")
  message("BURDEN STATS")
  message(paste0("Headline Estimate (Lambda=0.03): ", headline, " million"))
  message(paste0("Sensitivity Range (Lambda=0-0.064): ", lower, "–", upper, " million"))
  message("------------------------------------------------------")
  
} else {
  warning("Burden summary RDS not found. Run comparison script first.")
}
