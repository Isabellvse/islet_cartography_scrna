# Description -------------------------------------------------------------
# Here I plot results from differential gene expression of neighborhoods

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load_data ---------------------------------------------------------------
t2d <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/t2d_vs_nd.csv"))
pre <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/pre_vs_nd.csv"))

beta_up <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/t2d_vs_nd_deg_beta_down_vs_other.csv"))


# Beta --------------------------------------------------------------------
# (LAPTM5)Insulin granules undergo lysosomal degradation understress (stress-induced nascent granule degradation) https://www.nature.com/articles/s41467-019-11170-4
# (FCER1G) triggers allergic reaciton
# PTPRC (CD45) insulinitis https://www.researchgate.net/figure/CD45-staining-of-pancreas-A-Typical-example-of-Langerhans-islet-of-a-30-days-old-NOD_fig3_370302605
beta_up |> 
  dplyr::filter(padj <= 0.05 & log2FoldChange > 0) |> 
  dplyr::top_n(n = 20, wt = log2FoldChange)


df <- beta_up

# Geneset enrichment ------------------------------------------------------
# maybe also try go-term
# Get genesets to test
#https://biostatsquid.com/pathway-enrichment-analysis-tutorial-clusterprofiler/#step1
gene_sets_df <- msigdbr::msigdbr(species = 'Homo sapiens', collection = "C5")
term2gene <- gene_sets_df |> 
  dplyr::select(term = gs_name, gene = gene_symbol)

ranks <- df$log2FoldChange
names(ranks) <- df$gene_symbol
ranks_order <- sort(ranks, decreasing = T)

cp <- clusterProfiler::GSEA(geneList = ranks_order,
                      TERM2GENE = term2gene,
                      pvalueCutoff = 1, eps = 0)
test <- cp@result |> 
  dplyr::filter(p.adjust < 0.05)
