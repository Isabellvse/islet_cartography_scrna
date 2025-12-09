# Description -------------------------------------------------------------
# In this scripts, we will try to find a final annotation
# based on different metrics: reference mapping, module score, and gene expression 

# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
# Cell names map
markers_map <- vroom::vroom(here::here("islet_cartography_scrna/data/marker_database/cell_type_map.csv"))

# cell names
markers <- vroom::vroom(here::here("islet_cartography_scrna/data/marker_database/cell_type_map.csv")) |> 
  dplyr::mutate(n_char = nchar(new_cell_name)) |> 
  dplyr::arrange(desc(n_char)) |> 
  dplyr::pull("new_cell_name") |> 
  base::unique() |>
  base::paste(collapse = "|")


# Reference mapping
azimuth_annotation <- vroom::vroom(
  here::here("islet_cartography_scrna/data/annotate/files/azimuth_annotation.csv")
) |> 
  dplyr::mutate(
    # only keep annotations with high confidense 
    old_cell_name = dplyr::case_when(
      final_level_confidence < 0.6 ~ "unknown",
      .default = as.character(final_level_labels)
    ),
    old_cell_name = base::tolower(old_cell_name),
    # replace spaces, +, /, and all dash characters with _
    old_cell_name = stringr::str_replace_all(old_cell_name, "[ +/\\p{Pd}]", "_"),
    # remove parentheses
    old_cell_name = stringr::str_replace_all(old_cell_name, "[()]", ""),
    # collapse multiple underscores
    old_cell_name = stringr::str_replace_all(old_cell_name, "_+", "_")
  ) |> 
  dplyr::left_join(y = markers_map, by = "old_cell_name") |> 
  dplyr::select(barcode, azimuth_ref_anno = new_cell_name)

# Module scores
marker_gene_scores <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/marker_gene_scores.csv"))

# Manual annotation
manual_anno <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/manual_annotation.csv")) |> 
  dplyr::rename(old_cell_name = manual_annotation) |> 
  dplyr::left_join(y = markers_map, by = "old_cell_name") |> 
  dplyr::select(barcode, manual_anno = new_cell_name)

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
# **Since tf markers do  not cover all endocrine cells, I will not include them in the annotation**
module_long <- marker_gene_scores |> 
  tidyr::pivot_longer(cols = c(-barcode), names_to = "module", values_to = "score") |> 
  dplyr::left_join(module_columns, by = "module") |> 
  dplyr::filter(score > 0)

# Split module long into separate data frames for each database, and for each barcode, find the max module score
max_score <- module_long |> 
  (\(df) base::split(df, factor(df$database)))() |> 
  purrr::map(~ .x |> 
               dplyr::group_by(barcode) |> 
               dplyr::summarise(max = max(score))) |> 
  purrr::list_rbind(names_to = "database")

# Annotations
module_anno <- module_long |> 
  dplyr::left_join(max_score, by = c("barcode", "database")) |> 
  dplyr::mutate(is_cell_type = dplyr::case_when(score == max ~ 1,
                                                score != max ~ 0, 
                                                .default = 0)) |> 
  # Remove barcodes with no annotation
  dplyr::filter(is_cell_type == 1) |> 
  # Add cell type to those with an annotation
  dplyr::mutate(cell_type_anno = dplyr::case_when(is_cell_type == 1 ~ as.character(cell_type),
                                                  .default = as.character("unknown"))) |> 
  (\(df) base::split(df, factor(df$database)))() |> 
  purrr::imap( ~ .x  |> dplyr::select(barcode, cell_type_anno) |> 
                 # collapse annotations that score equally high (data.table for faster processing)
                 data.table::setDT() |> 
                 (\(df) df[, .(cell_type_anno = base::paste(cell_type_anno, collapse = "|")), by = barcode])() |> 
                 # rename columns with celltype annotation
                 dplyr::rename(!!.y := cell_type_anno) |> 
                 dplyr::rename_with(.cols = -c(barcode), function(x){paste0(x, "_anno")})) |> 
  # Join all dataframes
  purrr::reduce(dplyr::full_join) |> 
  # replace na values with unknown
  dplyr::mutate(dplyr::across(c(-barcode), ~.x |> tidyr::replace_na("unknown")))


test_4 <- module_anno[["panglao"]] |> 
  dplyr::group_by(barcode) |> 
  dplyr::summarise(Product = paste(cell_type_anno, collapse = "|"))
# some cell score max in more than 1 cell type
test <- module_anno[["panglao"]] |> 
  tidyr::pivot_wider(id_cols = c(barcode, database), 
                     names_from = cell_type_anno, 
                     values_from = is_cell_type) |> 
  dplyr::mutate(
  sum = base::rowSums(dplyr::pick(where(is.numeric)), na.rm = TRUE))

  # Join all dataframes
  purrr::reduce(dplyr::left_join) |> 
  # replace na values with unknown
  dplyr::mutate(dplyr::across(c(-barcode), ~.x |> tidyr::replace_na("unknown"))) |> 
  dplyr::rename_with(.cols = -c(barcode), function(x){paste0(x, "_anno")})

# Combine all annotations -------------------------------------------------
anno_comb <- purrr::reduce(list(manual_anno, module_anno), dplyr::full_join)



