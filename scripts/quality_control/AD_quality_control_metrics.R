# Description -------------------------------------------------------------
# Generate quality control metrics for each study, and save as a .csv file

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
## star_quality control ----
star_quality <- qs2::qs_read(here::here("islet_cartography_scrna/data/quality_control/star_quality_raw.qs2")) |> 
  # Add whether the method is droplet or plate based, ensure library prep is lower case
  purrr::modify_depth(1, ~ dplyr::mutate(., library_prep = base::tolower(library_prep), 
                                         platform = dplyr::case_when(library_prep %in% droplet_based ~ "droplet",
                                                                     library_prep %in% plate_based ~ "plate",
                                                                     library_prep %in% plate_based_bc ~ "plate_barcode")) %>% 
                        dplyr::relocate(platform, .after = "library_prep")) 

## Genes for qc ----
genes <- qs2::qs_read(here::here("islet_cartography_scrna/data/quality_control/mito_ribo_protein_genes.qs2"))
 
## Not empty barcodes
paths <- base::list.files(path = here::here("islet_cartography_scrna/data/quality_control/first_pass/barcode_ranks"), 
                          pattern = "*_not_empty", 
                          full.names = TRUE, 
                          recursive = TRUE) |> 
  purrr::set_names(~ stringr::str_extract(.x, pattern = "[^/]+_ic_\\d+_\\d+_\\d+"))

not_empty <- purrr::map(paths, ~vroom::vroom(., delim = ",", col_names = TRUE))

## mtx files ----
paths <- base::paste0("/work/scRNAseq/", base::names(star_quality), "/Preprocessed/*/Solo.out")
base::names(paths) <- purrr::map_chr(paths, ~stringr::str_extract(string = .x, pattern = "(?<=scRNAseq/)[^/]+"))

# quality control metrics -------------------------------------------------

purrr::imap(paths, function(path, name){
  base::message("Generating QC metrics for: ", name)
  quality_metrics_per_sample(path = path,
                          study_metadata = star_quality[[name]],
                          genes = genes,
                          study_name = name,
                          not_empty = not_empty)
  base::message("Finished ! generating QC metrics for: ", name)
  })

