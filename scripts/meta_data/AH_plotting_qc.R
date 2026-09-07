# Description -------------------------------------------------------------
# Plotting quality stats for datasets
# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
# Current meta data
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/obs.csv")) 

# raw meta data - was created when preparing anndata objects
meta_raw <- vroom::vroom(here::here("islet_cartography_scrna/data/metadata_harmonized/metadata_combined_harmonized.csv")) %>% 
  dplyr::rename(ic_id_dataset = ic_id_study)


# Combined plot of QC -----------------------------------------------------
# Dataset order
dataset_order <- meta |> 
  dplyr::distinct(ic_id_dataset) |> 
  dplyr::mutate(rank = as.numeric(stringr::str_remove(ic_id_dataset, "ic_"))) |> 
  dplyr::arrange(rank) |> 
  dplyr::pull(ic_id_dataset)

meta <- meta |> 
  dplyr::mutate(ic_id_dataset = factor(ic_id_dataset, levels = dataset_order))

p_ridges <- meta |> 
  dplyr::select(barcode, ic_id_dataset, n_umi, n_count, n_feature, 
                mitochondrial_fraction, coding_fraction, contrast_fraction) |> 
  tidyr::pivot_longer(cols = tidyselect::where(is.numeric)) |> 
  dplyr::mutate(name = factor(name, levels = rev(c("contrast_fraction", "coding_fraction", "mitochondrial_fraction",
                                                   "n_feature",
                                                   "n_count",
                                                   "n_umi")))) |> 
  ggplot2::ggplot(ggplot2::aes(x = value, y = ic_id_dataset)) +
  ggridges::geom_density_ridges(scale = 1, fill = "black") +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  ggplot2::facet_wrap(~name, scales = "free_x", nrow = 1) +
  my_theme() +
  ggplot2::theme(legend.position = "none")


n_cells_combined <- dplyr::full_join(
  dplyr::count(meta_raw, ic_id_dataset, name = "n_raw"),
  dplyr::count(meta, ic_id_dataset, name = "n_remaining"),
  by = "ic_id_dataset") |> 
  dplyr::transmute(ic_id_dataset,
                   Remaining = n_remaining,
                   Removed = n_raw - n_remaining) |> 
  tidyr::pivot_longer(-ic_id_dataset, names_to = "status", values_to = "n_cells") |> 
  dplyr::mutate(status = factor(status, levels = c("Removed", "Remaining")))

p_ncells <- ggplot2::ggplot(n_cells_combined, ggplot2::aes(n_cells, ic_id_dataset, fill = status)) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = c(Removed = "grey75", Remaining = "black")) +
  ggplot2::scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  my_theme() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank())

pdf(
  file = paste0(here::here("islet_cartography_scrna/data/paper_figures/plot/ridge_quality_control_plots.pdf")),
  height = 3,
  width = 5)
p_ridges + p_ncells + patchwork::plot_layout(widths = c(6, 1))
dev.off()


# cell type composition per dataset ---------------------------------------
all_cell_type_colors <- c(cell_type_colors, cell_type_broad_colors)
all_cell_type_colors <- all_cell_type_colors[!duplicated(names(all_cell_type_colors))]

pdf(
  file = paste0(here::here("islet_cartography_scrna/data/paper_figures/plot/cell_type_composition.pdf")),
  height = 3,
  width = 5)
meta |> 
  dplyr::select(barcode, ic_id_dataset, cell_type, cell_type_broad) |> 
  tidyr::pivot_longer(cols = c(cell_type, cell_type_broad),
                      names_to = "resolution", values_to = "cell_type_value") |> 
  dplyr::count(ic_id_dataset, resolution, cell_type_value) |> 
  dplyr::mutate(resolution = factor(resolution, levels = c("cell_type_broad", "cell_type"))) |> 
  ggplot2::ggplot(ggplot2::aes(x = n, y = ic_id_dataset, fill = cell_type_value)) +
  ggplot2::geom_col(position = "fill") +
  ggplot2::scale_fill_manual(values = all_cell_type_colors) +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5),
                              expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::facet_wrap(~resolution, nrow = 1) +
  my_theme() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                 legend.title = ggplot2::element_blank())
dev.off()
