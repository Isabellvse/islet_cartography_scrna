# Description -------------------------------------------------------------
# Add star alignment qc with ic ids

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/quality_control/"))
set.seed(1000)
# Load --------------------------------------------------------------------
## ic ids list ----
icid_list <- qs2::qs_read(here::here("islet_cartography_scrna/data/metadata/id_list.qs2"))

## Star quality control ----
files <- base::list.files(path = here::here("islet_cartography_scrna/data/star_quality"), pattern = ".csv$", full.names = T)
base::names(files) <- purrr::map_chr(files, ~ stringr::str_extract(string = .x, pattern = "(?<=alignment_qc_)[^/]+(?=\\.csv$)"))
star_quality <- purrr::imap(files, load_star_quality)

# Preprocess --------------------------------------------------------------
## Combine icids and star quality ----
star_quality_order <- star_quality[base::match(base::names(icid_list), base::names(star_quality))]
base::all.equal(base::names(icid_list), base::names(star_quality_order))
# Left join because we only want to keep samples that has been aligned.
star_id <- purrr::map2(icid_list, star_quality_order, dplyr::left_join)

# save --------------------------------------------------------------------
qs2::qs_save(star_id, here::here("islet_cartography_scrna/data/quality_control/star_quality_raw.qs2"))

