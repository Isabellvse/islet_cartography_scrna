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


# Overlap with harmonized labels ------------------------------------------
# Broad cell type annotaiton
meta |> 
  dplyr::mutate(agreement = as.numeric(study_cell_annotation_harmonized == cell_type_broad)) |> 
  dplyr::group_by(cell_type_broad) |>
  dplyr::summarise(total = dplyr::n(),
                   consistent = sum(agreement),
                   inconsistent = total-consistent,
                   pct_con = consistent / total,
                   pct_incon = inconsistent / total)

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
  file = paste0(here::here("islet_cartography_scrna/data/annotate/plot/annotation_agreement.pdf")),
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
