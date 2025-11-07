# Description -------------------------------------------------------------
# Get ribosomal, mitochondrial and protein-coding genes from gtf file

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/quality_control/"))
set.seed(1000)


# Load --------------------------------------------------------------------
# Only run this once
# gtf_data <- read_gtf_parallel("/work/islet_cartography_scrna/genome_files/gencode.v35.annotation.gtf", max_cores = 60, chunk_size = 20, max_size = 4 * 1024^3)

# Otherwise load this:
gtf_data <- qs2::qs_read(here::here("islet_cartography_scrna/genome_files/gencode.v35.annotation_gtf.qs2"))

# get genes ---------------------------------------------------------------
genes <- base::list(mito_genes = gtf_data |> dplyr::filter(seqname == "chrM") |> dplyr::pull(gene_id) |> unique(),
                    ribo_genes = gtf_data |> dplyr::filter(gene_type == "rRNA") |> dplyr::pull(gene_id) |> unique(),
                    protein_genes = gtf_data |> dplyr::filter(gene_type == "protein_coding") |> dplyr::pull(gene_id) |> unique())

# save --------------------------------------------------------------------
#qs2::qs_save(gtf_data, here::here("islet_cartography_scrna/genome_files/gencode.v35.annotation_gtf.qs2"))
qs2::qs_save(genes, here::here("islet_cartography_scrna/data/quality_control/mito_ribo_protein_genes.qs2"))


# Ensembl to gene map -----------------------------------------------------
gtf_data <- qs2::qs_read(here::here("islet_cartography_scrna/genome_files/gencode.v35.annotation_gtf.qs2"))

ensembl_gene <- gtf_data |> 
  dplyr::select(gene_id, gene_name) |> 
  dplyr::distinct() %>% 
  dplyr::mutate(gene_name_unique = make.unique(gene_name, sep = "_")) %>% 
  dplyr::select(-gene_name)

vroom::vroom_write(ensembl_gene, 
                   here::here("islet_cartography_scrna/genome_files/gene_id_map.csv"),
                   delim = ",", 
                   col_names = TRUE)


# Gene ensembl to entrez id map -------------------------------------------
entrez_gene <- ensembl_gene

entrez_gene$entrez_id <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                      keys = entrez_gene $gene_name_unique,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")

vroom::vroom_write(entrez_gene, 
                   here::here("islet_cartography_scrna/genome_files/gene_entrez_map.csv"),
                   delim = ",", 
                   col_names = TRUE)
# marker genes ------------------------------------------------------------
# Raw marker gene list
marker_dict <- list(
  cycling = c("UBE2C", "TOP2A", "CDK1", "BIRC5", "PBK", "CDKN3", "MKI67", "CDC20", "CCNB2", "CDCA3"),
  immune = c("ACP5", "APOE", "HLA-DRA", "TYROBP", "LAPTM5", "SDS", "FCER1G", "C1QC", "C1QB", "SRGN"),
  quiescent_stellate = c("RGS5", "C11orf96", "FABP4", "CSRP2", "IL24", "ADIRF", "NDUFA4L2", "GPX3", "IGFBP4", "ESAM"),
  endothelial = c("PLVAP", "RGCC", "ENG", "PECAM1", "ESM1", "SERPINE1", "CLDN5", "STC1", "MMP1", "GNG11"),
  schwann = c("NGFR", "CDH19", "UCN2", "SOX10", "S100A1", "PLP1", "TSPAN11", "WNT16", "SOX2", "TFAP2A"),
  activated_stellate = c("COL1A1", "COL1A2", "COL6A3", "COL3A1", "TIMP3", "TIMP1", "CTHRC1", "SFRP2", "BGN", "LUM"),
  epsilon = c("BHMT", "VSTM2L", "PHGR1", "TM4SF5", "ANXA13", "ASGR1", "DEFB1", "GHRL", "COL22A1", "OLFML3"),
  gamma = c("PPY", "AQP3", "MEIS2", "ID2", "GPC5-AS1", "CARTPT", "PRSS23", "ETV1", "PPY2", "TUBB2A"),
  delta = c("SST", "RBP4", "SERPINA1", "RGS2", "PCSK1", "SEC11C", "HHEX", "LEPR", "MDK", "LY6H"),
  ductal = c("SPP1", "MMP7", "IGFBP7", "KRT7", "ANXA4", "SERPINA1", "LCN2", "CFTR", "KRT19", "SERPING1"),
  acinar = c("REG1A", "PRSS1", "CTRB2", "CTRB1", "REG1B", "CELA3A", "PRSS2", "REG3A", "CPA1", "CLPS"),
  beta = c("IAPP", "INS", "DLK1", "INS-IGF2", "G6PC2", "HADH", "ADCYAP1", "GSN", "NPTX2", "C12orf75"),
  alpha = c("GCG", "TTR", "PPP1R1A", "CRYBA2", "TM4SF4", "MAFB", "GC", "GPX3", "PCSK2", "PEMT")
)

# Get all unique genes
all_genes <- sort(unique(unlist(marker_dict)))

# Initialize DataFrame with 0s
binary_markers <- as.data.frame(matrix(0, nrow = length(all_genes), ncol = length(marker_dict)))
rownames(binary_markers) <- all_genes
colnames(binary_markers) <- names(marker_dict)

# Fill with 1s where gene is a marker
for (celltype in names(marker_dict)) {
  binary_markers[marker_dict[[celltype]], celltype] <- 1
}

# Convery gene symbol to ensembl
binary_markers <- binary_markers |> 
  tibble::rownames_to_column("gene_name") |> 
  dplyr::left_join(y = ensembl_gene) |> 
  tidyr::drop_na() |> 
  dplyr::relocate(gene_id, .after = gene_name)

vroom::vroom_write(binary_markers, 
                   here::here("islet_cartography_scrna/data/binary_markers_gene_symbol.csv"),
                   delim = ",", 
                   col_names = TRUE)

vroom::vroom_write(binary_markers |> dplyr::select(-gene_name), 
                   here::here("islet_cartography_scrna/data/binary_markers_ensembl.csv"),
                   delim = ",", 
                   col_names = TRUE)
