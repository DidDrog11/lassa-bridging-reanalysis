## store your packages here

## core packages
library(renv)
library(dotenv)
library(conflicted) # make sure we know where function calls are coming from and which is preferred. See .Rprofile file for default preferences

## Project-Specific Packages
library(spOccupancy) # JSDM and Spatial Modelling
library(sf)
library(terra) # Efficient raster handling
library(rgbif) # Access GBIF API
library(coda) # MCMC

# Data Wrangling and Visualisation
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(lubridate)
library(tidyterra)

# Mapping and Environmental Data
library(geodata)
#install.packages("RcppArmadillo") # Might be needed for Luna
# install.packages('luna', repos='https://rspatial.r-universe.dev') # Not on CRAN
library(luna)
library(chirps)

# Epidemiology and Statistics
library(MASS) # Used for quasi-binomial GLMs or other stats
library(here) # File path management

## packages necessary for startup script
library(usethis)

## recommended
# library(fnmate) # quickly create functions


