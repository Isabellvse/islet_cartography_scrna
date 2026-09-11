# Description -------------------------------------------------------------
# Identify differentially expressed genes between disease groups for each cell type 

# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
base::source(here::here("islet_cartography_scrna/scripts/misc/marker_genes_functions.R"))
set.seed(1000)
vik <- khroma::color("vik")

library(furrr)
plan(multisession, workers = 120)

base_path <- here::here("islet_cartography_scrna/data/differential_genes_across_disease")
dir.create(base_path, showWarnings = FALSE)
dir.create(paste0(base_path, "/", "nd_t2d"), showWarnings = FALSE)
dir.create(paste0(base_path, "/", "nd_t2d/files"), showWarnings = FALSE)
dir.create(paste0(base_path, "/", "nd_t2d/plots"), showWarnings = FALSE)

# Load --------------------------------------------------------------------
df_paths <- base::list.files(path = here::here("islet_cartography_scrna/data/differential_genes_across_disease/nd_t2d"),
                             pattern = ".csv", 
                             full.names = T) |> 
  purrr::set_names(\(vec) (base::basename(vec) |>  
                             stringr::str_remove(".csv")))

df_paths |> 
  purrr::iwalk(\(x, idx) purrr::possibly(process_and_save)(path = x, name = idx, save_path = paste0(base_path, "/nd_t2d/files")))

## REMEMBER TO CHECK NUMBER OF DATASETS UNDERLYING THIS ANALYSIS 