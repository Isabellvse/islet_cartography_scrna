# Install helper packages if they are missing
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")

library(readr)
library(pak)
library(here)

# Install biocondictor version that fits R 4.4.3
install.packages("BiocManager")
BiocManager::install(version = "3.20")

# Load r packages
pkgs <- read_csv("/work/islet_cartography_scrna/practical_info/r_library.csv")

# Set library path
.libPaths(here::here("islet_cartography_scrna/r_library_2/"))

# Cache available package lists from cran and bioconductor
cran_pkgs <- available.packages(repos = "https://cran.r-project.org")[, "Package"]
bioc_pkgs <- BiocManager::available()

# Helper function to classify package origin (CRAN or bioconductor)
pkg_origin <- function(pkg) {
  if (pkg %in% cran_pkgs) {
    return("CRAN")
  } else if (pkg %in% bioc_pkgs) {
    return("Bioconductor")
  } else {
    return("Unknown")
  }
}

# Loop through packages and install accordingly
for (i in seq_len(nrow(pkgs))) {
  pkg <- pkgs$Package[i]
  ver <- pkgs$Version[i]
  origin <- pkg_origin(pkg)
  print(origin)
  if (origin == "Bioconductor") {
    message("Installing Bioconductor package: ", pkg)
    pak::pkg_install(paste0("bioc::", pkg))
  } else if (origin == "CRAN") {
    message("Installing CRAN package: ", pkg, " version ", ver)
    pak::pkg_install(paste0("cran/", pkg, "@", ver))
  } else {
    message("Skipping ", pkg, ": not found in CRAN or Bioconductor")
  }
}

.libPaths(here::here("islet_cartography_scrna/r_library/"))
write.csv(installed.packages()[, c("Package", "Version", "LibPath")],
          here("islet_cartography_scrna/practical_info/installed_packages_new.csv"), row.names = FALSE)

# Load old and new package lists
old <- read.csv(here("islet_cartography_scrna/practical_info/r_library.csv")) |>  dplyr::select(Package, Version, LibPath)
new <- read.csv(here("islet_cartography_scrna/practical_info/installed_packages_new.csv"))

# Find differences
diff <- merge(old, new, by = c("Package", "LibPath"), all = TRUE, suffixes = c("_old", "_new"))

# Show packages where versions differ
subset(diff, Version_old != Version_new)


# Downgrade packages
# Downgrade CRAN packages to old versions
pak::pkg_install("cran/circlize@0.4.16")
pak::pkg_install("cran/duckdb@1.4.1")
pak::pkg_install("cran/future@1.67.0")
pak::pkg_install("cran/future.apply@1.20.0")
pak::pkg_install("cran/ggbeeswarm@0.7.2")
pak::pkg_install("cran/ggplot2@4.0.0")
pak::pkg_install("cran/GlobalOptions@0.1.2")
pak::pkg_install("cran/gplots@3.2.0")
pak::pkg_install("cran/lme4@1.1-37")
pak::pkg_install("cran/parallelly@1.45.1")
pak::pkg_install("cran/progressr@0.17.0")
pak::pkg_install("cran/qs2@0.1.5")
pak::pkg_install("cran/RcppArmadillo@15.0.2-2")
pak::pkg_install("cran/reshape2@1.4.4")
pak::pkg_install("cran/reticulate@1.44.0")
pak::pkg_install("cran/RSQLite@2.4.3")
pak::pkg_install("cran/S7@0.2.0")
pak::pkg_install("cran/spatstat.explore@3.5-3")
pak::pkg_install("cran/spatstat.geom@3.6-0")
pak::pkg_install("cran/spatstat.random@3.4-2")
pak::pkg_install("cran/spatstat.univar@3.1-4")
pak::pkg_install("cran/uwot@0.2.3")
pak::pkg_install("cran/vroom@1.6.6")
pak::pkg_install("cran/XML@3.99-0.19")

# Downgrade Bioconductor packages (release-based, not version-pinned)
pak::pkg_install("bioc::GetoptLong")     # resolves to 1.0.5 in Bioc 3.20
pak::pkg_install("bioc::SeuratObject")   # resolves to 5.2.0 in Bioc 3.20
