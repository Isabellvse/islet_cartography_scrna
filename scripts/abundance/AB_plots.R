# Description -------------------------------------------------------------
# Here I plot differential abundance results

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load_data ---------------------------------------------------------------
t2d <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/t2d_vs_nd.csv"))
pre <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/pre_vs_nd.csv"))
df <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/de_analysis_overall/t2d_for_de_overall.csv"))


# Boxplot -----------------------------------------------------------------
pt2d = t2d |> 
  ggplot2::ggplot(ggplot2::aes(x = logFC, y = forcats::fct_reorder(nhood_annotation, logFC))) +
  ggrastr::geom_jitter_rast(aes(      
    colour = SpatialFDR <= 0.1), size = 0.1, shape = 20, raster.dpi = 1200) +
  ggplot2::geom_boxplot(outlier.shape = NA, fill = "transparent", color = "black", linewidth = 0.2) +
  ggplot2::labs(x = "Log2(T2D/ND)", y = "Annotated neighborhoods", colour = "SpatialFDR <= 10%") +
  ggplot2::scale_color_manual(values = c("FALSE" = "grey",
                                         "TRUE" = "#D9583B")) +
  ggplot2::geom_vline(xintercept = 0) +
  my_theme() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank())

ppre <- pre |> 
  ggplot2::ggplot(ggplot2::aes(x = logFC, y = forcats::fct_reorder(nhood_annotation, logFC))) +
  ggrastr::geom_jitter_rast(aes(      
    colour = SpatialFDR <= 0.1), size = 0.1, shape = 20, raster.dpi = 1200) +
  ggplot2::geom_boxplot(outlier.shape = NA, fill = "transparent", color = "black", linewidth = 0.2) +
  ggplot2::labs(x = "Log2(PRE/ND)", y = "Annotated neighborhoods", colour = "SpatialFDR <= 10%") +
  ggplot2::scale_color_manual(values = c("FALSE" = "grey",
                                         "TRUE" = "#E8B43F")) +
  ggplot2::geom_vline(xintercept = 0) +
  my_theme() 

pdf(
  file = here::here("islet_cartography_scrna/data/milo/plots/logfc_boxplot_t2d_pre.pdf"),
  height = 2,
  width = 5
)
(ppre + pt2d) + plot_layout(guides = 'collect')
dev.off()



# Vulcano plot ------------------------------------------------------------
pt2d = t2d |> 
  ggplot2::ggplot(ggplot2::aes(x = logFC, y = -log10(SpatialFDR))) +
  ggrastr::geom_point_rast(aes(      
    colour = SpatialFDR <= 0.1), size = 0.1, shape = 20, raster.dpi = 1200) +
  ggplot2::facet_wrap(~nhood_annotation, nrow = 1, scales = "free_y") +
  ggplot2::labs(x = "Log2(T2D/ND)", y = "-log10(SpatialFDR)", colour = "SpatialFDR <= 10%") +
  ggplot2::scale_color_manual(values = c("FALSE" = "grey",
                                         "TRUE" = "#D9583B")) +
  ggplot2::geom_vline(xintercept = 0) +
  my_theme()

ppre = pre |> 
  ggplot2::ggplot(ggplot2::aes(x = logFC, y = -log10(SpatialFDR))) +
  ggrastr::geom_point_rast(aes(      
    colour = SpatialFDR <= 0.1), size = 0.1, shape = 20, raster.dpi = 1200) +
  ggplot2::facet_wrap(~nhood_annotation, nrow = 1, scales = "free_y") +
  ggplot2::labs(x = "Log2(PRE/ND)", y = "-log10(SpatialFDR)", colour = "SpatialFDR <= 10%") +
  ggplot2::scale_color_manual(values = c("FALSE" = "grey",
                                         "TRUE" = "#E8B43F")) +
  ggplot2::geom_vline(xintercept = 0) +
  my_theme()


pdf(
  file = here::here("islet_cartography_scrna/data/milo/plots/volcano_t2d_pre.pdf"),
  height = 3,
  width = 15
)
(ppre / pt2d) + plot_layout(guides = 'collect')
dev.off()


# Beeswarmplot ------------------------------------------------------------
pt2d <- t2d |> 
  ggplot2::ggplot(ggplot2::aes(x = nhood_annotation, y = logFC)) +
  ggbeeswarm::geom_quasirandom(data=t2d[t2d$SpatialFDR > 0.1,], 
                               alpha=1, colour='grey50', size = 0.1) +
  ggbeeswarm::geom_quasirandom(data=t2d[t2d$SpatialFDR <= 0.1,], 
                               aes(colour=logFC), size = 0.1) +
  ggplot2::scale_color_distiller(palette = "RdBu", direction = -1, limits = c(-2, 2)) +
  ggplot2::coord_flip() +
  ggplot2::labs(x ="Annotated neighborhoods", y = "Log2(T2D/ND)", colour = "Log2FC") +
  my_theme() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank())

ppre <- pre |> 
  ggplot2::ggplot(aes(x = nhood_annotation, y = logFC)) +
  ggbeeswarm::geom_quasirandom(data=pre[pre$SpatialFDR > 0.1,], 
                               alpha=1, colour='grey50', size = 0.1) +
  ggbeeswarm::geom_quasirandom(data=pre[pre$SpatialFDR <= 0.1,], 
                               aes(colour=logFC), size = 0.1) +
  ggplot2::scale_color_distiller(palette = "RdBu", direction = -1, limits = c(-2, 2)) +
  ggplot2::coord_flip() +
  ggplot2::labs(x ="Annotated neighborhoods", y = "Log2(PRE/ND)", colour = "Log2FC") +
  my_theme() 

pdf(
  file = here::here("islet_cartography_scrna/data/milo/plots/beeswarm_t2d_pre.pdf"),
  height = 2,
  width = 5
)
(ppre + pt2d) + plot_layout(guides = 'collect')
dev.off()


# Number of cells in each neighborhood ------------------------------------
table <- df |> 
  dplyr::group_by(nhood_id) |> 
  dplyr::summarise(n_cells = dplyr::n()) |> 
  dplyr::summarise(median = round(median(n_cells), 2),
                   mean = round(mean(n_cells), 2),
                   min = round(min(n_cells), 2),
                   max = round(max(n_cells), 2))
pcell <- df |> 
  dplyr::group_by(nhood_id) |> 
  dplyr::summarise(n_cells = dplyr::n()) |> 
  ggplot2::ggplot(ggplot2::aes(x = n_cells)) +
  ggplot2::geom_histogram() +
  ggplot2::annotate(geom = 'table',
                    x = 7000, 
                    y = 100,
                    label=list(table),
                    table.theme = gridExtra::ttheme_minimal(
                      base_size = 4)) +
  ggplot2::labs(title = "Number of cells",
                x = "Total number of cells\nper neighborhood",
                y = "Count") +
  my_theme()


# Number of samples -------------------------------------------------------
table <- df |> 
  dplyr::select(nhood_id, ic_id_platform_adjusted_sample) |> 
  dplyr::distinct() |> 
  dplyr::group_by(nhood_id) |> 
  dplyr::summarise(n_samples = dplyr::n()) |> 
  dplyr::summarise(median = round(median(n_samples), 2),
                   mean = round(mean(n_samples), 2),
                   min = round(min(n_samples), 2),
                   max = round(max(n_samples), 2))


psample <- df |> 
  dplyr::select(nhood_id, ic_id_platform_adjusted_sample) |> 
  dplyr::distinct() |> 
  dplyr::group_by(nhood_id) |> 
  dplyr::summarise(n_samples = dplyr::n()) |> 
  ggplot2::ggplot(ggplot2::aes(x = n_samples)) +
  ggplot2::geom_histogram() +
  ggplot2::labs(title = "Number of samples",
                x = "Total number of samples\nper neighborhood",
                y = "Count") +
  my_theme()

# Purity ------------------------------------------------------------------
table <- df |> 
  dplyr::select(nhood_id, nhood_annotation_frac) |> 
  dplyr::distinct() |> 
  dplyr::summarise(median = round(median(nhood_annotation_frac), 3),
                   mean = round(mean(nhood_annotation_frac), 3),
                   min = round(min(nhood_annotation_frac), 3),
                   max = round(max(nhood_annotation_frac), 3))

ppure <- df |> 
  dplyr::select(nhood_id, nhood_annotation_frac) |> 
  dplyr::distinct() |> 
  ggplot2::ggplot(ggplot2::aes(x = nhood_annotation_frac)) +
  ggplot2::geom_histogram(bins = 100) +
  ggplot2::annotate(geom = 'table',
                    x = 0.2, 
                    y = 300,
                    label=list(table),
                    table.theme = gridExtra::ttheme_minimal(
                      base_size = 4)) +
  ggplot2::labs(title = "Purity",
                x = "Fraction of cell types\nin neighborhood",
                y = "Count") +
  my_theme()
