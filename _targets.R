## For more information about targets see https://books.ropensci.org/targets
## Load your packages, e.g. library(targets).
source("./packages.R")
# load environment variables
dotenv::load_dot_env()

# Set target options with targets::tar_option_set() - useful for tuning performance and setting up pipelines on HPCs

## Load functions stored in /R 
targets::tar_source()

## logical collection of targets that does some task
example_tar_plan <- tar_plan(
  
  # tar_target(targets_are_nouns, functions_are_verbs(arg)) 
  
)

example_tar_plan2 <- tar_plan(
  
   tar_target(data, get_data(arg)),
   tar_target(data_cleaned, clean_data(data))
  
)

# add tar_plan outputs to this list to have targets make them.
list(
  example_tar_plan,
  example_tar_plan2
)
