## store your packages here

## core packages
library(renv)
library(dotenv)
library(conflicted) # make sure we know where function calls are coming from and which is preferred. See .Rprofile file for default preferences

## Project-Specific Packages
library(spOccupancy) # JSDM and Spatial Modelling
library(sf)
library(terra) # Efficient raster handling

# Data Wrangling and Visualisation
library(dplyr)
library(tidyr)
library(ggplot2)

# Epidemiology and Statistics
library(MASS) # Used for quasi-binomial GLMs or other stats
library(here) # File path management

## packages necessary for startup script
library(usethis)

## recommended
# library(fnmate) # quickly create functions


