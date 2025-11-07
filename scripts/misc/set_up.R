# Install 'here' globally if it's missing
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")   # goes to your default user library
}

# library path
.libPaths(here::here("islet_cartography_scrna/r_library/"))

# source files
base::source(here::here("islet_cartography_scrna/scripts/misc/package_load.R"))
base::source(here::here("islet_cartography_scrna/scripts/misc/functions.R"))
base::source(here::here("islet_cartography_scrna/scripts/misc/misc.R"))
