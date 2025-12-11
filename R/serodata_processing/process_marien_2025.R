library(tidyverse)
library(here)

# 1. Define the Guinea Data (Mariën et al. / Faranah)
# Study Period: Jan-March 2014
# Note: Yarawalia coordinates provided in DMS, others in DDM.

guinea_raw <- tibble(Village = c("Brissa", "Dalafinani", "Damania", "Sokourala", "Sonkonia", "Yarawalia"),
                     NumTested = c(234, 211, 203, 213, 235, 210),
                     NumPositive = c(165, 168, 181, 189, 216, 156),
                     # Brissa: 10 13.010 -> 10 + 13.010/60 = 10.21683
                     # Dalafinani: 10 08.590 -> 10 + 8.590/60 = 10.14317
                     # Damania: 09 48.410 -> 9 + 48.410/60 = 9.80683
                     # Sokourala: 10 03.407 -> 10 + 3.407/60 = 10.05678
                     # Sonkonia: 09 54.763 -> 9 + 54.763/60 = 9.91272
                     # Yarawalia: 9 57 0 -> 9 + 57/60 = 9.95000
                     Latitude = c(10.21683, 10.14317, 9.80683, 10.05678, 9.91272, 9.95000),
                     # Longitude (West is Negative)
                     # Brissa: 10 41.326 -> -(10 + 41.326/60) = -10.68877
                     # Dalafinani: 10 36.303 -> -(10 + 36.303/60) = -10.60505
                     # Damania: 10 51.796 -> -(10 + 51.796/60) = -10.86327
                     # Sokourala: 10 39.950 -> -(10 + 39.950/60) = -10.66583
                     # Sonkonia: 10 47.888 -> -(10 + 47.888/60) = -10.79813
                     # Yarawalia: 10 43 59.99 -> -(10 + 43/60 + 59.99/3600) = -10.73333
                     Longitude = c(-10.68877, -10.60505, -10.86327, -10.66583, -10.79813, -10.73333))

# 2. Format for Pipeline
guinea_clean <- guinea_raw |>
  mutate(Country = "Guinea",
         Admin1 = "Faranah",
         Year = 2014, 
         Source = "https://doi.org/10.1093/infdis/jiaf308",
         Human_Random_Survey = TRUE) |>
  dplyr::select(Town_Region = Admin1, Village, Year, Latitude, Longitude, Country, NumTestAb = NumTested, NumPosAb = NumPositive, Human_Random_Survey, Source)

csv_path <- here("data", "raw", "human_sero_larocs.csv")
write_csv(guinea_clean, csv_path)
