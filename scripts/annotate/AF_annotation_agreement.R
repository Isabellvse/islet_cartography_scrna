# Description -------------------------------------------------------------
# In this scripts, we will try to find a final annotation
# based on different metrics: reference mapping, module score, and gene expression 

# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
# Reference mapping
azimuth_annotation <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/azimuth_annotation.csv"))
# Module scores
marker_gene_scores <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/marker_gene_scores.csv"))
# Marker genes
markers <- vroom::vroom(here::here("islet_cartography_scrna/data/marker_database/cell_type_map.csv")) |> 
  dplyr::mutate(n_char = nchar(new_cell_name)) |> 
  dplyr::arrange(desc(n_char)) |> 
  dplyr::pull("new_cell_name") |> 
  base::unique() |>
  base::paste(collapse = "|")

# Manual annotation
manual_anno <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/manual_annotation.csv"))


# Preprocessing -----------------------------------------------------------
# Make database and celltype columns for marker genes
module_columns <- marker_gene_scores |> 
  tidyr::pivot_longer(c(-barcode), names_to = "module", values_to = "score") |> 
  dplyr::select(module) |> 
  dplyr::distinct() |> 
  dplyr::mutate(database = stringr::str_remove(module, "_UCell"),
                database = stringr::str_remove(database, markers),
                database = stringr::str_remove(database, "_$"),
                cell_type = stringr::str_extract(module, markers))

# Annotation based on module scores ---------------------------------------
# Make module score long, and add database and cell type
module_long <- marker_gene_scores |> 
  tidyr::pivot_longer(cols = c(-barcode), names_to = "module", values_to = "score") |> 
  dplyr::left_join(module_columns, by = "module") 

# Split module long into seperate dataframes for each database, and for each barcode, find the max module score

max_score <- module_long |> 
  (\(df) base::split(df, factor(df$database)))() |> 
  purrr::map(~ .x |> 
               dplyr::group_by(barcode) |> 
               dplyr::summarise(max = max(score))) |> 
  purrr::list_rbind(names_to = "database")

# add max module score to module_long, by database and barcode.
# if a score is equal to max score and higher than 0, the cell gets 1, if it is not it gets 0
module_anno <- module_long |> 
  dplyr::left_join(max_score, by = c("barcode", "database")) |> 
  dplyr::mutate(is_cell_type = dplyr::case_when(score == max & score > 0 ~ 1,
                                             score != max ~ 0, 
                                             .default = 0)) |> 
  # Remove barcodes with no annotation
  dplyr::filter(is_cell_type == 1) |> 
  # Add cell type to those with an annotation
  dplyr::mutate(cell_type_anno = dplyr::case_when(is_cell_type == 1 ~ as.character(cell_type),
                                                  .default = as.character("unknown"))) |> 
  (\(df) base::split(df, factor(df$database)))() |> 
  # rename columns with celltype annotation
  purrr::imap( ~ .x  |> dplyr::select(barcode, cell_type_anno) |>  
                 dplyr::rename(!!.y := cell_type_anno)) |> 
  # Join all dataframes
  purrr::reduce(dplyr::full_join) |> 
  # replace na values with unknown
  dplyr::mutate(dplyr::across(c(-barcode), ~.x |> tidyr::replace_na("unknown")))

  

