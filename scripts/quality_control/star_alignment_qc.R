# Description -------------------------------------------------------------
# Add star alignment qc with ic ids

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/quality_control/"))

# Load --------------------------------------------------------------------
# load ic ids and split into a list according to study name
icid_list <- qs::qread(here::here("islet_cartography_scrna/data/metadata/meta_combined_2.qs")) |>
  dplyr::select(ic_id, ic_id_study, ic_id_donor, ic_id_sample, name, sample, library_prep) |>
  (\(df) base::split(df, factor(df$name)))()

# load star quality csv
# Get paths
files <- base::list.files(path = here::here("islet_cartography_scrna/data/star_quality"), pattern = ".csv$", full.names = T)

# Add names
base::names(files) <- purrr::map_chr(files, ~ stringr::str_extract(string = .x, pattern = "(?<=alignment_qc_)[^/]+(?=\\.csv$)"))

# Load all csv files
df_list <- purrr::imap(files, load_star_quality)

# Make sure lists have the same order
df_list_order <- df_list[base::match(base::names(icid_list), base::names(df_list))]

# Check correct order
base::all.equal(base::names(df_list_order), base::names(icid_list))

# Add ic ids to star quality control --------------------------------------
# Keeping only samples that have been aligned
star_qc_list <- purrr::map2(df_list_order, icid_list, dplyr::left_join)
 
# Sanity check that we dont loose information about samples
base::all.equal(base::names(df_list_order), base::names(star_qc_list))
purrr::map2(df_list_order, star_qc_list, ~base::all.equal(base::nrow(.x), base::nrow(.y)))
