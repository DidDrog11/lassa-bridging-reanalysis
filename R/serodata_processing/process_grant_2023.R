install.packages("tidyverse")
library(tidyverse)
library(sf)
library(terra)
library(here)

# 1. Load GADM Data for Sierra Leone
# You downloaded this in get_data.R, assuming it's in data/raw/gadm/
gadm_path <- here("data", "raw", "gadm")
sl_gadm_files <- list.files(gadm_path, pattern = "SLE.*\\.rds", full.names = TRUE)
sl_adm1 <- readRDS(grep("_1_", sl_gadm_files, value = TRUE)[1]) 
sl_adm1_sf <- st_as_sf(sl_adm1)

# 2. Clean the PEER Data Table
peer_clean <- PEER_Locations |>
  # A. Fill Village names down to the empty 'W' rows
  fill(Village, `Community Size*`, Participants, `Seropositivity N (%)`, .direction = "down") |>
  # B. Identify Lat vs Lon rows
  mutate(Type = ifelse(grepl("N", `GPS Coordinates`), "Lat", "Lon"),
    # Extract just the number (removes 'N', 'W', spaces)
    Value = as.numeric(gsub("[^0-9.]", "", `GPS Coordinates`))) |>
  # C. Pivot to one row per village
  pivot_wider(id_cols = c(Village, `Community Size*`, Participants, `Seropositivity N (%)`),
              names_from = Type,
              values_from = Value) |>
  # D. Format Coordinates (West = Negative)
  # The raw numbers are likely positive (e.g. 11.36), so we force them negative.
  mutate(Lon = -abs(Lon)) |>
  # E. Extract Counts
  mutate(N_Positive = as.numeric(str_extract(`Seropositivity N (%)`, "^\\d+")),
         N_Tested = Participants) |>
  filter(!is.na(Lat) & !is.na(Lon)) |>
  dplyr::select(Village, Lat, Lon, N_Tested, N_Positive)

# 3. Reverse Geocode (Assign Admin Regions)
peer_sf <- st_as_sf(peer_clean, coords = c("Lon", "Lat"), crs = 4326)
peer_geocoded <- st_join(peer_sf, sl_adm1_sf)

# 4. Final Formatting
peer_final <- peer_geocoded |>
  st_drop_geometry() |>
  select(Village, NAME_1) |>
  left_join(peer_clean, by  = "Village") |>
  mutate(Country = "Sierra Leone",
         Admin1 = NAME_1, 
         Year = 2016, 
         Source = "https://doi.org/10.1371/journal.pntd.0010938",
         Human_Random_Survey = TRUE) |>
  dplyr::select(Town_Region = Admin1, Village, Year, Latitude = Lat, Longitude = Lon, Country, NumTestAb = N_Tested, NumPosAb = N_Positive, Human_Random_Survey, Source)

# 5. Save
write.csv(peer_final, here("data", "raw", "human_sero_PEER.csv"), row.names = FALSE)
