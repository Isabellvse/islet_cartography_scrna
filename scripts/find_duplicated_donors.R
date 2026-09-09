# Description -------------------------------------------------------------
# Find duplicated donors across datasets, prioritize keeping those made with 10x,
# or those with most cell (if none has been used for 10x)

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)
library(glmmTMB)
library(emmeans)

# Load --------------------------------------------------------------------
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/obs.csv"))

# Find duplicated donors --------------------------------------------------
meta |> 
  dplyr::mutate(
    library_prep_clean = dplyr::if_else(grepl("10x", library_prep, ignore.case = TRUE), "10x", "non-10x"),
    priority = dplyr::if_else(library_prep_clean == "10x", 1, 2)
  ) |> 
  dplyr::group_by(ic_id_donor_overall, ic_id_dataset, priority, library_prep_clean) |> 
  dplyr::tally() |> 
  dplyr::group_by(ic_id_donor_overall) |> 
  dplyr::filter(dplyr::n_distinct(ic_id_dataset) > 1) |> 
  dplyr::arrange(ic_id_donor_overall, priority, dplyr::desc(n)) |> 
  dplyr::slice(-1) |> # remove the "best" donors
  dplyr::ungroup() |> 
  dplyr::select(ic_id_donor_overall, ic_id_dataset) |> 
  vroom::vroom_write(here::here("islet_cartography_scrna/data/duplicated_donors.csv"))
