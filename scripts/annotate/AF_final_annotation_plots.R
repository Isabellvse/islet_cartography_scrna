# Description -------------------------------------------------------------
# Plotting results from manual annotation

# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)


# Misc --------------------------------------------------------------------
markers <- list(
  "acinar" = c("PRSS1", "PRSS2", "CPA1"),
  "acinar_reg_plus" = c("REG3A", "REG3G", "REG1B"),
  "alpha" = c("GCG", "TTR", "MAFB"),
  "beta" = c("INS", "IAPP", "INS-IGF2"),
  "cycling" = c("TOP2A","MKI67", "CCNB1"),
  "delta" = c("SST", "HHEX", "LY6H"),
  "ductal" = c("KRT19","CFTR","HNF1B"),
  "ductal_mucin" = c("MUC1", "TFF1", "TFF2"),
  "endmt" = c("COL3A1", "DDR2", "HIF1A"),
  "endothelial" = c("PECAM1", "PLVAP", "FLT1"),
  "endothelial_islet" = c("PASK", "ESM1", "LAMA4"),
  "epsilon" = c("GHRL", "PHGR1", "RBP4"),
  "gamma" = c("PPY", "ARX","ETV1"),
  "mast" = c("S100A4", "RGS10", "LTC4S"),
  "myeloid" = c("LYZ", "HLA-DRA", "CD68"),
  "schwann" = c("CRYAB", "S100B", "PMP22"),
  "stellate_a" = c("TIMP1", "COL1A2", "VCAN"),
  "stellate_q" = c("RGS5", "FABP4", "ADIRF"))

markers_df <- enframe(markers, name = "cell_type", value = "gene_symbol") |>
  unnest(gene_symbol)

vik <- khroma::color("vik")
# Load --------------------------------------------------------------------
deg <- list.files(
  path = here::here("islet_cartography_scrna/data/annotate/dg_onevsother"),
  pattern = "vs_allother.csv",
  recursive = TRUE,
  full.names = TRUE
) |>
  # keep only files NOT containing "pd_barcode"
  keep(~ !stringr::str_detect(.x, "pd_barcode")) |> 
  purrr::map(~vroom::vroom(.x) |> dplyr::filter(padj <= 0.05 & pct_expr_c1 > 0.5 & log2FoldChange >=1)) |> 
  dplyr::bind_rows()

mean_markers <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/mean_marker_expression.csv"))
var <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/var_marker_expression.csv"))
pct <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/genes_pct_expressed.csv"))

# Preprocess --------------------------------------------------------------
pct_df <- pct |> 
  tidyr::pivot_longer(tidyselect::where(is.numeric), names_to = "manual_annotation", values_to = "pct") |> 
  dplyr::mutate(manual_annotation = stringr::str_remove(manual_annotation, "pct_expr_cluster_"))

mean_df <- mean_markers |> 
  tidyr::pivot_longer(tidyselect::where(is.numeric), names_to = "gene_symbol", values_to = "gene_exp") |> 
  dplyr::mutate(gene_symbol = stringr::str_remove(gene_symbol, "\\...*")) |> 
  dplyr::left_join(y = markers_df,
                   relationship = "many-to-many") |> 
  dplyr::left_join(y = pct_df) |> 
  dplyr::group_by("manual_annotation") |> 
  dplyr::mutate(z_score = base::scale(gene_exp)) |> 
  dplyr::ungroup()

rng <- range(mean_df$z_score, na.rm = TRUE)

pdf(
  file = here::here("islet_cartography_scrna/data/annotate/plot/dotplot_markergenes.pdf"),
  height = 4,
  width = 8
)
mean_df |>  
  ggplot2::ggplot(ggplot2::aes(x = gene_symbol, y = manual_annotation)) +
  ggplot2::geom_point(ggplot2::aes(size = pct, fill = gene_exp), color = "black", shape = 21) +
  ggplot2::scale_size("Score", range = c(0, 4)) +
    ggplot2::scale_fill_gradientn(
      colors = vik(256),
      limits = c(-3, 3),
      values = scales::rescale(c(-3, 0, 3)),
      oob = scales::squish
    ) +
  ggplot2::facet_wrap(~cell_type, nrow = 1,  scales = "free_x", labeller = labeller(
    cell_type = function(x)
      ggplot2::label_wrap_gen(width = 10)(
        stringr::str_to_sentence(gsub("_", " ", x))
      ))) +
  ggplot2::scale_y_discrete(
    labels = function(y)(
        stringr::str_to_sentence(gsub("_", " ", y))
      )
  ) +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                 panel.spacing = unit(0.1, "lines"))
  
dev.off()
