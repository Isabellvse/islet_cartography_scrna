# Description -------------------------------------------------------------
# Here I plot results from differential gene expression of neighborhoods

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load_data ---------------------------------------------------------------
t2d <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/t2d_vs_nd.csv"))
pre <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/pre_vs_nd.csv"))

# df data
df <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/de_analysis_overall/t2d_for_de_overall.csv"))

# Load deg files 
dg <- Sys.glob(file.path(here::here("islet_cartography_scrna/data/milo/files"),"t2d_vs_nd_deg*_vs_other.csv")) |> 
  purrr::set_names(~stringr::str_remove(base::basename(.), "t2d_vs_nd_deg_") |> stringr::str_remove("_vs_other.csv"))
dg_list <- purrr::map(dg, vroom::vroom)

# Mean gene expression per neighborhood
mean_exp_nhood <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/mean_counts_per_nhood.csv"))

# Overall diagnostics -----------------------------------------------------
# Number of cells in each neighborhood
pfrac <- df |> 
  dplyr::group_by(nhood_id, disease_harmonized) |>
  dplyr::summarise(n = n()) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(nhood_id) |> 
  dplyr::mutate(sum = sum(n),
                frac = n/sum) |> 
  #tidyr::pivot_wider(id_cols = nhood_id, values_from = frac, names_from = disease_harmonized, values_fill = 0) |> 
  ggplot2::ggplot(ggplot2::aes(x = disease_harmonized, y = nhood_id, fill = frac)) +
  ggplot2::geom_tile() +
  my_theme()
  
pdf(
  file = here::here("islet_cartography_scrna/data/milo/plots/n_cell_neighborhood_condition.pdf"),
  height = 2,
  width = 2
)
pfrac
dev.off()


pfrac <- df |> 
  dplyr::group_by(nhood_id, ic_id_platform_adjusted_sample) |>
  dplyr::summarise(n = n()) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(nhood_id) |> 
  dplyr::mutate(sum = sum(n),
                frac = n/sum) |> 
  #tidyr::pivot_wider(id_cols = nhood_id, values_from = frac, names_from = ic_id_platform_adjusted_sample, values_fill = 0) |> 
  ggplot2::ggplot(ggplot2::aes(x = ic_id_platform_adjusted_sample, y = nhood_id, fill = frac)) +
  ggplot2::geom_tile() +
  my_theme()

df |>
  dplyr::distinct(nhood_id, ic_id_platform_adjusted_sample, disease_harmonized) |>
  dplyr::count(nhood_id, disease_harmonized) |>
  dplyr::group_by(nhood_id) |>
  dplyr::mutate(frac = n / sum(n)) |> 
  ggplot(aes(x = factor(nhood_id), y = frac, fill = disease_harmonized)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Preprocess data ---------------------------------------------------------
# Upregulated genes
dg_sig <- dg_list |> 
  purrr::map(~dplyr::filter(., padj <= 0.05 & log2FoldChange >= 1.5))

# Vulcano -----------------------------------------------------------------
p_vul <- dg_list |> 
  purrr::imap(\(df, group){
    output <- df |>   
      tidyr::drop_na() |> 
      ggplot2::ggplot(ggplot2::aes(x = log2FoldChange, y = -log10(padj))) +
      ggrastr::geom_point_rast(aes(      
        colour = padj <= 0.05), size = 0.1, shape = 20, raster.dpi = 1200) +
      ggplot2::labs(x = base::paste0("Log2(", group, "/all)"), y = "-log10(FDR)", colour = "FDR <= 5%") +
      ggplot2::scale_color_manual(values = c("FALSE" = "grey",
                                             "TRUE" = "#D9583B")) +
      ggplot2::geom_vline(xintercept = 0) +
      my_theme()
    return(output)
  }())

pdf(
  file = here::here("islet_cartography_scrna/data/milo/plots/volcano_deg_t2d_neighborhoods.pdf"),
  height = 2,
  width = 3
)
p_vul
dev.off()

t2d |> 
  ggplot2::ggplot(ggplot2::aes(x = logFC, y = -log10(SpatialFDR))) +
  ggrastr::geom_point_rast(aes(      
    colour = SpatialFDR <= 0.1), size = 0.1, shape = 20, raster.dpi = 1200) +
  ggplot2::facet_wrap(~nhood_annotation, nrow = 1, scales = "free_y") +
  ggplot2::labs(x = "Log2(T2D/ND)", y = "-log10(SpatialFDR)", colour = "SpatialFDR <= 10%") +
  ggplot2::scale_color_manual(values = c("FALSE" = "grey",
                                         "TRUE" = "#D9583B")) +
  ggplot2::geom_vline(xintercept = 0) +
  my_theme()




# Heatmap -----------------------------------------------------------------
# MAke a loop with heatmaps?
# See what cells are in these neighborhoods? 
# These should probably also be divided into nd and t2d?
# what is the distance from the neighbors?


# Heatmap with logFC og nhood ---------------------------------------------
# Order of neighborhoods
t2d_beta <- t2d |> 
  dplyr::filter(nhood_annotation == "beta")

ranks <- t2d_beta$logFC
names(ranks) <- t2d_beta$index_cell
ranks_order <- sort(ranks)

ranks_color <- rep(NA_character_, length(ranks_order))
names(ranks_color) <- names(ranks_order)
ranks_color[ranks_order > 0.5] <- "red"
ranks_color[ranks_order <= 0.5] <- "blue"
ranks_color[names(ranks_order) %in% (t2d_beta |> dplyr::filter(SpatialFDR > 0.1) |> dplyr::pull(index_cell))] <- "grey"

# Beta up
beta_up_sig <- dg_sig |> 
  purrr::pluck("beta_up")

# # z-score scale
# mat <- mean_exp_nhood |> 
#   dplyr::filter(gene_symbol %in% beta_up_sig$gene_symbol) |> 
#   tidyr::pivot_longer(-gene_symbol) |> 
#   tidyr::pivot_wider(id_cols = name, names_from=gene_symbol, values_from=value) |> 
#   dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ (. - mean(.)) / (sd(.)))) |> 
#   tidyr::pivot_longer(-name, names_to = "gene_symbol", values_to = "value") |> 
#   tidyr::pivot_wider(id_cols = gene_symbol, names_from=name, values_from=value) |> 
#   tibble::column_to_rownames("gene_symbol") |> 
#   as.matrix()

mat <- mean_exp_nhood |> 
  dplyr::filter(gene_symbol %in% beta_up_sig$gene_symbol) |> 
  tidyr::pivot_longer(-gene_symbol) |> 
  tidyr::pivot_wider(id_cols = name, names_from=gene_symbol, values_from=value) |> 
  dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ (. - min(.)) / (max(.)-min(.)))) |> 
  tidyr::pivot_longer(-name, names_to = "gene_symbol", values_to = "value") |> 
  tidyr::pivot_wider(id_cols = gene_symbol, names_from=name, values_from=value) |> 
  tibble::column_to_rownames("gene_symbol") |> 
  as.matrix()

# Reorder matrix by logFC from nhood
mat <- mat[, names(ranks_order)]

# Gene groups
gene_logFC <- beta_up_sig$log2FoldChange
names(gene_logFC) <- beta_up_sig$gene_symbol
gene_group <- ifelse(gene_logFC > 0, "Up", "Down")
# Order of genes should match matrix
gene_group <- gene_group[rownames(mat)]

# Point plot for neighborhoods
col_ha <- ComplexHeatmap::HeatmapAnnotation(
  logFC = ComplexHeatmap::anno_points(
    ranks_order,
    pch = 16,
    size = unit(2, "mm"), 
    gp = grid::gpar(col = ranks_color)
  )
)

breaks <- seq(0, 1, length.out = 100)
colors <- viridis::inferno(100)

ComplexHeatmap::Heatmap(
  mat,
  name = "Mean expression\n (z-score)",
  top_annotation = col_ha,
  show_row_dend = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  col = circlize::colorRamp2(breaks, colors),
  row_split = gene_group,
  show_row_names = FALSE,
  show_column_names = FALSE,
  border_gp = grid::gpar(col = "black"),
  column_title = "Neighborhoods",
  column_title_side = "bottom",
  use_raster = TRUE
)


# Geneset enrichment analysis ---------------------------------------------
gsea <- dg_list |> 
  purrr::map(\(df){
    # rank genes
    ranks <- df$log2FoldChange
    names(ranks) <- df$gene_symbol
    ranks_order <- sort(ranks, decreasing = T)
    
    # GSEA
    gse <- clusterProfiler::gseGO(geneList=ranks_order, 
                                  ont ="BP", 
                                  keyType = "SYMBOL", 
                                  pvalueCutoff = 0.05, 
                                  verbose = TRUE, 
                                  OrgDb = org.Hs.eg.db::org.Hs.eg.db)
    
    gse_simple <- clusterProfiler::simplify(gse)
    
    return(gse_simple)
  }())

# Filter data ----
gsea_sig <- gsea |> 
  purrr::map(\(res){
    output <- res@result |>
      dplyr::filter(!is.na(p.adjust)) |>
      dplyr::filter(p.adjust <= 0.05) |>
      dplyr::filter(abs(NES) >= 1.5) |>
      dplyr::filter(stringr::str_count(core_enrichment, "/") >= 1) |> #Filter out gene with only 1 core enrichment gene
      dplyr::arrange(dplyr::desc(NES))
    return(output)
  }())

# save --------------------------------------------------------------------
qs2::qs_save(gsea, here::here("islet_cartography_scrna/data/milo/objects/gsea_res.qs2"))
qs2::qs_save(gsea_sig, here::here("islet_cartography_scrna/data/milo/objects/gsea_res_sig.qs2"))


















