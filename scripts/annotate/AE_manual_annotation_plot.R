# Description -------------------------------------------------------------
# Plotting results from manual annotation (differential expressed genes)

# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

## list of endocrine markers ----
# https://pubmed.ncbi.nlm.nih.gov/18797460/
# https://www.nature.com/articles/s41419-021-03603-0#Sec2
# https://www.nature.com/articles/s41419-021-03920-4

endo_markers <- list(
  endothelial = c("PLVAP", "CLDN5", "PECAM1", "VWF", "CDH5", "KDR"),
  
  islet_endothelial_cells = c(
    "ACE", "PASK", "F2RL3", "ESM1", "CXCR4", "UNC5B", "LAMA4",
    "CREM", "COL13A1", "NKX2-3", "ANGPTL2", "THBS1"
  ),
  endmt = c("S100A4", "APOE", "TGFBI", "C3AR1", "LGALS3")
)

acinar_markers <- list(
  Acinar_i = c(
    "RBPJL", "CHRM3", "LRIG1", "INSR",
    "FOXP2", "CHN2", "DTNA", "SDK1",
    "MAP3K5", "CAMK1D"
  ),
  
  Acinar_REG_plus = c(
    "REG3A", "REG3G", "REG1B", "REG1A"
  ),
  
  Acinar_s = c(
    "CPB1", "CPA1", "PRSS3", "PRSS1",
    "AMY2A", "CELA2A", "CELA3A", "CELA3B",
    "CELA3A", "CELA3B", "CTRB1", "CTRB2",
    "CLPS", "PNLIP", "SPINK1", "CTRC",
    "CPA2", "ANXA4"
  )
)


## Define colors for heatmap ----
myCol <- colorRamp::colorRampPalette(c('#004B7A', 'white', '#A83708'))(100)
myBreaks <- base::seq(-1.5, 1.5, length.out = 100)
col_fun <- circlize::colorRamp2(myBreaks, myCol)

# Load --------------------------------------------------------------------
c_22 <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/deg_wald_22_subcluster_vs_23.csv"))
c_14 <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/deg_wald_14_subcluster_0_vs_1.csv"))
c_17 <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/deg_wald_17_subcluster_0_vs_1.csv"))
c_1718 <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/deg_wald_17_vs_18.csv"))

# Cluster 22 heatmap ------------------------------------------------------
# Generate dataframe with group of genes
gene_groups <- tibble::tibble(
  gene_symbol = base::unlist(endo_markers),
  group = base::rep(base::names(endo_markers), base::lengths(endo_markers))
)

# Create vector of genes of interest
en_mak <- endo_markers |> unlist() |> unname() |> unique()

# Matrix with log2 fold changes
mat <- c_22 |>
  dplyr::filter(gene_symbol %in% en_mak) |>
  dplyr::select(gene_symbol, comparison, log2FoldChange) |>
  tidyr::pivot_wider(id_cols = gene_symbol,
                     names_from = comparison,
                     values_from = log2FoldChange) |>
  tibble::column_to_rownames("gene_symbol") |>
  as.matrix()

# Order gene groups
gene_groups_ordered <- gene_groups[match(rownames(mat), gene_groups$gene_symbol), ]
row_split_vector <- gene_groups_ordered$group

# Significance star
sig <- c_22 |>
  dplyr::filter(gene_symbol %in% en_mak) |>
  dplyr::mutate(padj = dplyr::case_when(padj <= 0.01 ~ "*",
                                        .default = "")) |>
  dplyr::select(gene_symbol, comparison, padj) |>
  tidyr::pivot_wider(id_cols = gene_symbol,
                     names_from = comparison,
                     values_from = padj) |>
  tibble::column_to_rownames("gene_symbol") |>
  as.matrix()

pdf(here::here("islet_cartography_scrna/data/annotate/plot/cluster22_logfc_heatmap.pdf"), 
               width = 3, height = 6)
ComplexHeatmap::Heatmap(mat,
        col = col_fun,
        name = "Log2FC",
        row_split = row_split_vector, # Use the corrected, ordered vector
        cluster_rows = TRUE, # Hierarchical clustering *within* each split
        cluster_columns = FALSE,
        cluster_row_slices = FALSE, # This is good if you want to cluster within splits
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(sig[i, j], x, y, gp = grid::gpar(fontsize = 10))},
        rect_gp = grid::gpar(col = "black", lwd = 1)
)
dev.off()

# Cluster 14 heatmap ------------------------------------------------------
# Generate dataframe with group of genes
gene_groups <- tibble::tibble(
  gene_symbol = base::unlist(acinar_markers),
  group = base::rep(base::names(acinar_markers), base::lengths(acinar_markers))
)

# Create vector of genes of interest
ac_mak <- acinar_markers |> unlist() |> unname() |> unique()

# Matrix with log2 fold changes
mat <- c_14 |>
  dplyr::filter(gene_symbol %in% ac_mak) |>
  dplyr::select(gene_symbol, comparison, log2FoldChange) |>
  tidyr::pivot_wider(id_cols = gene_symbol,
                     names_from = comparison,
                     values_from = log2FoldChange) |>
  tibble::column_to_rownames("gene_symbol") |>
  as.matrix()

# Order gene groups
gene_groups_ordered <- gene_groups[match(rownames(mat), gene_groups$gene_symbol), ]
row_split_vector <- gene_groups_ordered$group

# Significance star
sig <- c_14 |>
  dplyr::filter(gene_symbol %in% ac_mak) |>
  dplyr::mutate(padj = dplyr::case_when(padj <= 0.01 ~ "*",
                                        .default = "")) |>
  dplyr::select(gene_symbol, comparison, padj) |>
  tidyr::pivot_wider(id_cols = gene_symbol,
                     names_from = comparison,
                     values_from = padj) |>
  tibble::column_to_rownames("gene_symbol") |>
  as.matrix()

pdf(here::here("islet_cartography_scrna/data/annotate/plot/cluster14_logfc_heatmap.pdf"), 
    width = 3, height = 6)
ComplexHeatmap::Heatmap(mat,
                        col = col_fun,
                        name = "Log2FC",
                        row_split = row_split_vector, # Use the corrected, ordered vector
                        cluster_rows = TRUE, # Hierarchical clustering *within* each split
                        cluster_columns = FALSE,
                        cluster_row_slices = FALSE, # This is good if you want to cluster within splits
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid::grid.text(sig[i, j], x, y, gp = grid::gpar(fontsize = 10))},
                        rect_gp = grid::gpar(col = "black", lwd = 1)
)
dev.off()

# cluster 17 --------------------------------------------------------------
c_17 |>
  dplyr::filter(padj <= 0.05)

# cluster 17 and 18 -------------------------------------------------------
# https://www.gastrojournal.org/article/S0016-5085(20)35399-3/fulltext?referrer=https%3A%2F%2Fpubmed.ncbi.nlm.nih.gov%2F
muc_genes <- c("MUC5B", "TFF1", "TFF2", "TFF3", "CRISP3", "CFTR", "SLC4A4", "SCTR")

test <- c_1718 |>
  dplyr::filter(gene_symbol %in% muc_genes) %>% 
  dplyr::mutate(sig = dplyr::case_when(padj <= 0.01 ~ "*",
                                        .default = "")) %>% 
  dplyr::relocate(sig)
  dplyr::filter(padj <= 0.01 & abs(log2FoldChange) >= 1) |>
  dplyr::arrange(desc(log2FoldChange), padj, desc(pct_expr_c1))

# The smaller subtype (accounting for 1% of the total ductal cells) 
# was characterized by higher expression levels of genes linked to mucous 
# secretion such as MUC5B (hereafter, MUC5B+ ductal cells); 
# the trefoil factor genes TFF1, TFF2, and TFF3; 
# and the cysteine rich secretory protein 3 CRISP3 (Figure 2C and D). 
# The other ductal subtype, by contrast, showed higher expression levels of 
# classical ductal markers such CFTR, SLC4A4, and SCTR (Figure 2C and D).3 