# Description -------------------------------------------------------------
# Here I plot differential abundance results

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load_data ---------------------------------------------------------------
t2d <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/t2d_vs_nd.csv"))
pre <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/pre_vs_nd.csv"))


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
  ggplot2::ggplot(aes(x = nhood_annotation, y = logFC)) +
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
