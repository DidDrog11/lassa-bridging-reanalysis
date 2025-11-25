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
  arha_rds_url <- "https://raw.githubusercontent.com/DidDrog11/arenavirus_hantavirus/main/data/database/Project_ArHa_database_2025-09-18.rds"
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

west_africa_region <- gadm(country = west_africa_countries$iso3c,  level = 1, path = data_dir)

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

# 4. Acquire Epidemiological Data (Tasks C & D) ---------------------------------

## C. Human Seroprevalence Data
# This is likely contained within the Lassa_Occurrence_Data file from the original paper.
# e.g. seroprevalence_raw <- read_excel(here("data", "raw", "Lassa_Occurrence_Data.xlsx"), 
#                                      sheet = "Cleaned_Lassa_Literature")
# # Save as a clean CSV
# write.csv(seroprevalence_raw, here(data_dir, "raw_human_seroprevalence.csv"), row.names = FALSE)


## D. Reported Case Data (Validation Data - CRITICAL NEW DATA)
# Define sources for national/sub-national reported Lassa Fever cases.
# This often requires manual search and cleaning of published reports/WHO data.
# Placeholder for the data input script/API call.
# e.g. who_case_reports <- read.csv(here("data", "raw", "who_lf_case_reports.csv"))


# 5. Review and Finalize -------------------------------------------------------

message("Data acquisition script complete.")

# END OF SCRIPT