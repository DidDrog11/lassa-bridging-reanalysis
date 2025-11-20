# load packages
source(here::here("packages.R"))

## update the readme file 
usethis::edit_file(here::here("README.md"))

## update .env file
usethis::edit_file(here::here(".env"))
dotenv::load_dot_env()

# update license - usethis provides many options for licenses
usethis::use_mit_license(copyright_holder = Sys.getenv("copyright_holder"))

## setup renv for package management
renv::init()

# use renv::hydrate to get any packages already installed on your machine
#renv::hydrate()
# use renv::status to see if any packages need to be installed
#renv::status()
# use renv::install to install missing packages
#renv::install("openssl")

