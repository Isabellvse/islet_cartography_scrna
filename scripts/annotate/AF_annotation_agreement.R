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
  dplyr::pull("new_cell_name") |> 
  base::unique() |> 
  base::paste(sep = "_",collapse = "|")
# Manual annotation
manual_anno <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/manual_annotation.csv"))


# Annotation based on module scores ---------------------------------------
module <- marker_gene_scores |> 
  tidyr::pivot_longer(cols = c(-barcode), names_to = "module", 
                      values_to = "score") |> 
  dplyr::mutate(database = stringr::str_remove(module, pattern = markers),
                database = stringr::str_remove(module, pattern = "__UCell"),
                cell_type = stringr::str_extract(module, pattern = markers))

