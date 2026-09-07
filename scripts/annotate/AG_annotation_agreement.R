# Description -------------------------------------------------------------
# cell type agreement

# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
base::source(here::here("islet_cartography_scrna/scripts/misc/marker_genes_functions.R"))
set.seed(1000)

iridescent <- khroma::color("iridescent")

# Load --------------------------------------------------------------------
# Current meta data
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/obs.csv")) 

# Broad annotation --------------------------------------------------------
confusion <- meta |>
  dplyr::count(cell_type_broad, 
               study_cell_annotation_harmonized, 
               name = "n") |>
  tidyr::complete(
    cell_type_broad,
    study_cell_annotation_harmonized,
    fill = list(n = 0)
  ) |> 
  dplyr::group_by(cell_type_broad) |>
  dplyr::mutate(pct = n / sum(n)) |>
  dplyr::ungroup() 

# of the cells I labeled X, what fraction were called Y by the study
pdf(
  file = paste0(here::here("islet_cartography_scrna/data/annotate/plot/annotation_agreement_cell_type_broad.pdf")),
  height = 2,
  width = 3)
confusion |>
  ggplot2::ggplot(ggplot2::aes(x = study_cell_annotation_harmonized, y = cell_type_broad, fill = pct)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradientn(colors = iridescent(5), limits=c(0,1)) +
  ggplot2::labs(
    x = "Study annotation",
    y = "Cell type broad",
    fill = "Fraction of\n cell type broad"
  ) +
  ggplot2::scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_title(x)) +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()


# cell type ---------------------------------------------------------------
confusion <- meta |>
  dplyr::count(cell_type, 
               study_cell_annotation_harmonized, 
               name = "n") |>
  tidyr::complete(
    cell_type,
    study_cell_annotation_harmonized,
    fill = list(n = 0)
  ) |> 
  dplyr::group_by(cell_type) |>
  dplyr::mutate(pct = n / sum(n)) |>
  dplyr::ungroup() 

# of the cells I labeled X, what fraction were called Y by the study
pdf(
  file = paste0(here::here("islet_cartography_scrna/data/annotate/plot/annotation_agreement_cell_type.pdf")),
  height = 2,
  width = 3)
confusion |>
  ggplot2::ggplot(ggplot2::aes(x = study_cell_annotation_harmonized, y = cell_type, fill = pct)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradientn(colors = iridescent(5), limits=c(0,1)) +
  ggplot2::labs(
    x = "Study annotation",
    y = "Cell type",
    fill = "Fraction of\n cell type"
  ) +
  ggplot2::scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_title(x)) +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()


# Per dataset Broad annotation --------------------------------------------------------
confusion <- meta |>
  dplyr::count(ic_id_dataset,
               cell_type_broad, 
               study_cell_annotation_harmonized, 
               name = "n") |>
  tidyr::complete(
    ic_id_dataset,
    cell_type_broad,
    study_cell_annotation_harmonized,
    fill = list(n = 0)
  ) |> 
  dplyr::group_by(ic_id_dataset, cell_type_broad) |>
  dplyr::mutate(pct = n / sum(n)) |>
  dplyr::ungroup()

pdf(
  file = paste0(here::here("islet_cartography_scrna/data/annotate/plot/annotation_agreement_cell_type_broad_dataset.pdf")),
  height = 5,
  width = 6)
confusion |>
  ggplot2::ggplot(ggplot2::aes(x = study_cell_annotation_harmonized, y = cell_type_broad, fill = pct)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradientn(colors = iridescent(5), limits = c(0, 1)) +
  ggplot2::facet_wrap(~ ic_id_dataset) +
  ggplot2::labs(
    x = "Study annotation",
    y = "Cell type broad",
    fill = "Fraction of\nmy label"
  ) +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

# Per dataset annotation --------------------------------------------------------
confusion <- meta |>
  dplyr::count(ic_id_dataset,
               cell_type, 
               study_cell_annotation_harmonized, 
               name = "n") |>
  tidyr::complete(
    ic_id_dataset,
    cell_type,
    study_cell_annotation_harmonized,
    fill = list(n = 0)
  ) |> 
  dplyr::group_by(ic_id_dataset, cell_type) |>
  dplyr::mutate(pct = n / sum(n)) |>
  dplyr::ungroup()

pdf(
  file = paste0(here::here("islet_cartography_scrna/data/annotate/plot/annotation_agreement_cell_type_dataset.pdf")),
  height = 6,
  width = 6)
confusion |>
  ggplot2::ggplot(ggplot2::aes(x = study_cell_annotation_harmonized, y = cell_type, fill = pct)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradientn(colors = iridescent(5), limits = c(0, 1)) +
  ggplot2::facet_wrap(~ ic_id_dataset) +
  ggplot2::labs(
    x = "Study annotation",
    y = "Cell type",
    fill = "Fraction of\nmy label"
  ) +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

# % of cells that are mislabeled -------------------------------------------
meta |> 
  dplyr::mutate(inconsistent = as.numeric(cell_type_broad != study_cell_annotation_harmonized)) |> 
  dplyr::group_by(cell_type_broad) |> 
  dplyr::summarise(total = dplyr::n(),
                   sum = sum(inconsistent),
                   pct = sum / total) |> 
  ggplot2::ggplot(ggplot2::aes(x = pct, y = cell_type_broad, fill = cell_type_broad)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = cell_type_broad_colors) +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_title(x)) +
  my_theme() +
  ggplot2::theme(legend.position = "none")

meta |> 
  dplyr::mutate(inconsistent = as.numeric(cell_type != study_cell_annotation_harmonized)) |> 
  dplyr::group_by(cell_type) |> 
  dplyr::summarise(total = dplyr::n(),
                   sum = sum(inconsistent),
                   pct = sum / total) |> 
  ggplot2::ggplot(ggplot2::aes(x = pct, y = cell_type, fill = cell_type)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = cell_type_colors) +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_title(x)) +
  my_theme() +
  ggplot2::theme(legend.position = "none")
