################################################################################
## SCRIPT: get_data.R
##
## PURPOSE: Acquire all necessary raw data inputs for the Lassa Spillover 
##          Reanalysis. This includes historical rodent records, new GBIF
##          occurrences for competitor species, environmental rasters, 
##          and validation data (seroprevalence and case reports).
##
## OUTPUTS: Saves all raw data files to the ./data/raw/ subdirectory.
################################################################################

# 1. Setup and Package Loading -------------------------------------------------

# Load core packages (defined in packages.R)
source(here::here("packages.R"))

# Define output directory structure
data_dir <- here("data", "raw")
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}
proc_dir <- here("data", "processed")
if (!dir.exists(proc_dir)) {
  dir.create(proc_dir, recursive = TRUE)
}
# 2. Acquire Rodent Occurrence Data (Task A) ------------------------------------

## A. Original Paper's Data (Manual Input/Local Read)
# Read in the original Mastomys natalensis presence data
if(!file.exists(here(data_dir, "mastomys_presences_original.csv"))) {
  mastomys_presences_original <- read_tsv("https://raw.githubusercontent.com/54481andrew/pathogen-spillover-forecast/refs/heads/master/Reservoir_Layer/Data/Mastomys_natalensis_presences.csv")
  write_csv(mastomys_presences_original, here(data_dir, "mastomys_presences_original.csv"))
  lasv_sequences_original <- read_tsv("https://raw.githubusercontent.com/54481andrew/pathogen-spillover-forecast/refs/heads/master/Reservoir_Layer/Data/Rodents_Genbank_Oct_2020.csv")
  write_csv(lasv_sequences_original, here(data_dir, "lasv_sequences_original.csv"))
} else {
  mastomys_presences_original <- read_csv(here(data_dir, "mastomys_presences_original.csv"))
  lasv_sequences_original <- read_csv(here(data_dir, "lasv_sequences_original.csv"))
}

## B. Data Source ArHa
if(!file.exists(here(data_dir, "Project_ArHa_database.rds"))) {
  arha_rds_url <- "https://raw.githubusercontent.com/DidDrog11/arenavirus_hantavirus/main/data/database/Project_ArHa_database_2025-12-02.rds"
  local_file_path <- here("data", "raw", "Project_ArHa_database.rds")
  download.file(arha_rds_url, local_file_path, mode = "wb") 
  arha_db <- read_rds(local_file_path)
} else {
  arha_db <- read_rds(here(data_dir, "Project_ArHa_database.rds"))
}

## C. GBIF Data
species_to_obtain <- arha_db$pathogen |> 
  filter(str_detect(pathogen_species_cleaned, "lassa")) |> 
  filter(number_positive >= 1) |> 
  distinct(host_record_id) |> 
  left_join(arha_db$host, by  = "host_record_id") |> 
  drop_na(host_species) |> 
  distinct(host_species, gbif_id) |> 
  filter(!str_detect(host_species, "Bandicota|Rattus argentiventer|Rattus exulans")) # Remove rodents from outside WA

west_africa_countries <- countrycode::codelist |>
  filter(str_detect(region23, "Western Africa")) |>
  select(country = country.name.en, iso3c, iso2c) |>
  filter(!str_detect(iso3c, "CPV|SHN")) # Remove island nations

west_africa_region <- gadm(country = west_africa_countries$iso3c, level = 1, path = data_dir)
writeVector(west_africa_region, here(data_dir, "west_africa_shapefile.rds"))
west_africa_ext <- ext(west_africa_region)

# gbif_job <- occ_download(# Filters for quality and scope
#   pred_in("hasCoordinate", TRUE),
#   pred_gte("year", 1960),
#   pred_lte("year", 2025),
#   pred_in("country", west_africa_countries$iso3c),
#   # Filter by all taxon keys
#   pred_in("taxonKey", species_to_obtain$gbif_id),
#   format = "SIMPLE_CSV"
# )
# 
# gibf_data <- occ_download_wait('0000676-251120083545085') %>%
#   occ_download_import()
# GBIF Citation
# GBIF.org (20 November 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.k3pv43
gbif_data <- read_tsv(here(data_dir, "gbif_data.csv"))

# D. Explicit urban sampling
urban_lit_path <- here(data_dir, "rodent_urban_literature.csv")

# E. West African Rodents dataset
if(!file.exists(here(data_dir, "WA_rodents_database.rds"))) {
  wa_rodents_rds_url <- "https://raw.githubusercontent.com/DidDrog11/data_for_gbif/main/data_raw/rodent_data.rds"
  local_file_path <- here("data", "raw", "WA_rodents_database.rds")
  download.file(wa_rodents_rds_url, local_file_path, mode = "wb") 
  wa_rodent_db <- read_rds(local_file_path)
} else {
  wa_rodent_db <- read_rds(here(data_dir, "WA_rodents_database.rds"))
}

# 3. Acquire Environmental Predictors (Task B) ----------------------------------
start <- "2001-01-01"
end <- "2025-06-30" # Use the end of June as the last data point for environmental conditions

# Bio1 = Annual Mean Temperature (Proxy for Tmu)
# Bio12 = Annual Precipitation (Proxy for Pmu)
wc_data <- worldclim_global(var = "bio",  res = 0.5, path = data_dir)
# Select the two layers of interest
Tmu_proxy_wc <- wc_data[["wc2.1_30s_bio_1"]]
Pmu_proxy_wc <- wc_data[["wc2.1_30s_bio_12"]]

# NDVI data
ndvi_product_mod <- "MOD13C2" # Terra/MODIS Monthly NDVI
ndvi_product_myd <- "MYD13C2" # Aqua/MODIS Monthly NDVI

username_eosdis <- Sys.getenv("EOSDIS_USER")
password_eosdis <- Sys.getenv("EOSDIS_PWD")

# mf_mod <- luna::getNASA(product = ndvi_product_mod, start = start, end = end, aoi = west_africa_ext,
#                         download = TRUE, overwrite = FALSE, path = here(data_dir, "ndvi"), username = username_eosdis, password = password_eosdis)
# mf_myd <- luna::getNASA(product = ndvi_product_myd, start = start, end = end, aoi = west_africa_ext,
#                         download = TRUE, overwrite = FALSE, path = here(data_dir, "ndvi"),, username = username_eosdis, password = password_eosdis)

# Precipitation data
chirps_dir <- here("data", "raw", "chirps")
if (!dir.exists(chirps_dir)) {
  dir.create(chirps_dir, recursive = TRUE)
}

month_starts <- seq(ymd(start), ymd(end), by = "month")
month_ends <- month_starts + months(1) - days(1)
month_ends[length(month_ends)] <- end

monthly_precip_list <- list()

for (i in 1:length(month_starts)) {
  start_m <- month_starts[i]
  end_m <- month_ends[i]
  dates_range <- c(format(start_m, "%Y-%m-%d"), format(end_m, "%Y-%m-%d"))
  
  output_filename <- here(chirps_dir, paste0("chirps_daily_", format(start_m, "%Y%m"), ".tif"))
  
  if (!file.exists(output_filename)) {
    message(paste("Downloading daily data for:", format(start_m, "%Y-%m")))
    tryCatch({
      # Request daily data for the current month.
      # object = SpatExtent: Returns SpatRaster cropped to the extent.
      # resolution = 0.05
      daily_stack <- get_chirps(object = vect(west_africa_ext, crs = "EPSG:4326"), dates = dates_range, server = "CHC", resolution = 0.05, as.raster = TRUE)
      
      # Calculate the MONTHLY SUM (total precipitation)
      monthly_sum <- app(daily_stack, fun = "sum")
      names(monthly_sum) <- paste0("Psum_", format(start_m, "%Y%m"))
      
      writeRaster(monthly_sum, output_filename, overwrite = TRUE)
    }, error = function(e) {
      warning(paste("Error downloading CHIRPS for", format(start_m, "%Y%m"), ":", e$message))
    })
  } else {
    message(paste("Skipping:", output_filename, "already exists."))
  }
}

# Landcover data
lc_dir <- here("data", "raw", "landcover")
if (!dir.exists(lc_dir)) {
  dir.create(lc_dir, recursive = TRUE)
}

landcover_product <- "MCD12Q1"

years <- year(ymd(start)):year(ymd(end))

# for (year in years) {
#   start_date <- paste0(year, "-01-01")
#   end_date <- paste0(year, "-12-31")
#   
#   mf_lc <- luna::getNASA(product = landcover_product, start = start_date, end = end_date, aoi = west_africa_ext, 
#                          download = TRUE, overwrite = FALSE, path = lc_dir, 
#                          username = username_eosdis, password = password_eosdis)
# }

# Elevation
# Global raster
elev_rast <- elevation_global(res = 0.5, path = data_dir)

# Population density
pop_dir <- here("data", "raw", "pop")
# Manually downloaded from https://hub.worldpop.org/geodata/listing?id=77
# Unconstrained individual countries 2000-2020 UN adjusted (only 2020 used)

# Nighttime Lights (VIIRS)
ntl_dir <- here(data_dir, "ntl")
if (!dir.exists(ntl_dir)) dir.create(ntl_dir)
# Manually 

# 4. Acquire Epidemiological Data (Tasks C & D) ---------------------------------

## A. Human Seroprevalence Data
human_sero_url <- "https://raw.githubusercontent.com/54481andrew/pathogen-spillover-forecast/master/Pathogen_Layer/Data/Cleaned_Lassa_Literature.csv"
human_sero_file <- here(data_dir, "human_seroprevalence_basinski.csv")

if (!file.exists(human_sero_file)) {
  # Note: The original file is tab-separated despite the .csv extension
  human_sero_raw <- readr::read_tsv(human_sero_url)
  
  write_csv(human_sero_raw, human_sero_file)
} else {
  human_sero_raw <- readr::read_tsv(human_sero_file)
}

all_additional <- readr::read_delim(here("data", "raw", "human_seroprevalence_updated.csv"), 
                                    delim = ";", 
                                    show_col_types = FALSE)

## B. Reported Case Data (Validation Data)
# 1. Previously compiled data
# Define sources for national/sub-national reported Lassa Fever cases.
lasv_cases_moore <- "https://raw.githubusercontent.com/mooresea/lassa-model/main/data/case_reports_Lassa_CASES_full.csv"
lassa_cases_file_moore <-  here(data_dir, "lassa_cases_moore.csv")

if (!file.exists(lassa_cases_file_moore)) {
  lassa_cases_raw <- readr::read_csv(lasv_cases_moore)
  
  write_csv(lassa_cases_raw, lassa_cases_file_moore)
} else {
  lassa_cases_raw <- readr::read_tsv(lassa_cases_file_moore)
}

# 2. Annual binary ADM2 cases
ncdc_files <- list.files(here(data_dir), pattern = "ncdc_affected_lgas_20.*\\.csv", full.names = TRUE)

if (length(ncdc_files) > 0) {
  # 1. Bind them all together
  ncdc_combined <- purrr::map_df(ncdc_files, readr::read_csv, show_col_types = FALSE)
  # 3. Aggregate: Identify "Ever-Reported" LGAs
  # If an LGA appears in ANY year, it is Endemic (Status = 1)
  # We group by Code/State/LGA to ensure uniqueness
  ncdc_status_positive <- ncdc_combined |> 
    group_by(LGA_Code) |> 
    summarise(State = first(State),
              LGA = first(LGA),
              Years_Reported = n(),
              Last_Reported = max(Year),
              Status = 1,
              .groups = "drop")
  
  # 4. Merge with Full GADM List to get the "0"s (Non-Endemic)
  # Load Nigeria GADM (Admin 2)
  nga_gadm_file <- list.files(here(data_dir, "gadm"), pattern = "NGA.*_2_.*\\.rds", full.names = TRUE)[1]
  
  if (!is.na(nga_gadm_file)) {
    nga_polys <- readRDS(nga_gadm_file)
    
    all_lgas <- data.frame(LGA_Code = nga_polys$GID_2,
                           State = nga_polys$NAME_1,
                           LGA = nga_polys$NAME_2) |>
      left_join(ncdc_status_positive, by = "LGA_Code") |>
      mutate(Status = ifelse(is.na(Status), 0, 1),
             State = coalesce(State.x, State.y),
             LGA = coalesce(LGA.x, LGA.y)) |>
      dplyr::select(LGA_Code, State, LGA, Status, Years_Reported, Last_Reported)
    
    write.csv(all_lgas, here(proc_dir, "validation_binary_lga.csv"), row.names = FALSE)
    message(paste("Binary Validation Set: ", sum(all_lgas$Status), "Endemic LGAs out of", nrow(all_lgas)))
    
  } else {
    warning("Nigeria GADM file not found. Cannot generate binary zeros.")
  }
} else {
  warning("No NCDC annual CSVs found.")
}


# 3. NCDC ADM1 Case Counts (State-Level Check) ---
# Reads the 'ncdc_lassa_data_adm1' file
adm1_path <- here(data_dir, "ncdc_lassa_data_adm1.csv")

if (file.exists(adm1_path)) {
  # Read (handling potential messy headers if manual)
  ncdc_adm1 <- readr::read_csv(adm1_path, show_col_types = FALSE) |>
    # Filter out totals rows
    filter(states != "Total") |>
    # Aggregate Confirmed Cases by State
    group_by(states) |>
    summarise(Total_Confirmed = sum(confirmed, na.rm = TRUE), .groups="drop") |>
    rename(State = states)
  
  write.csv(ncdc_adm1, here(proc_dir, "validation_cases_nigeria_adm1.csv"), row.names = FALSE)
}


# 5. Basinski rasters for comparison --------------------------------------
basinski_host_layer <- "https://raw.githubusercontent.com/54481andrew/pathogen-spillover-forecast/master/Reservoir_Layer/Figures_Fits/reservoir_v7/Mn_pa_nboots_25_nbg_Same_tc_1_mllr_2_lmt_7/Reservoir_Layer_Mn_pa_nboots_25_nbg_Same_tc_1_mllr_2_lmt_7.tif"
basinski_pathogen_layer <- "https://raw.githubusercontent.com/54481andrew/pathogen-spillover-forecast/master/Pathogen_Layer/Figures_Fits/pathogen_v7/pa_nboots_25_1_mllr_3_lmt_7_ambi_NA_mintest_5/Lassa_Layer_pa_nboots_25_1_mllr_3_lmt_7_ambi_NA_mintest_5.tif"
basinski_combined_risk_layer <- "https://raw.githubusercontent.com/54481andrew/pathogen-spillover-forecast/master/Human_LASV_Incidence/Figures_Fits/human_v7_ambiNA_mint_5/Combined_Risk_Layer.tif"

basinski_path <- here("data", "raw", "basinski_comparison")

if(!file.exists(here(basinski_path, "host_raster.tif"))) {
  
  basinski_host_raster <- rast(basinski_host_layer)
  writeRaster(basinski_host_raster, here(basinski_path, "host_raster.tif"))
  basinski_pathogen_raster <- rast(basinski_pathogen_layer)
  writeRaster(basinski_pathogen_raster, here(basinski_path, "pathogen_raster.tif"))
  basinski_combined_raster <- rast(basinski_combined_risk_layer)
  writeRaster(basinski_combined_raster, here(basinski_path, "lasv_risk_raster.tif"))
}

# 6. Review and Finalize -------------------------------------------------------

message("Data acquisition script complete.")

# END OF SCRIPT