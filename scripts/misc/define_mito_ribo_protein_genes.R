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

# feature type ------------------------------------------------------------
gene_type <- gtf_data |> 
  dplyr::select(gene_symbol = gene_name, feature_type = gene_type) |> 
  dplyr::distinct()

vroom::vroom_write(gene_type, 
                   here::here("islet_cartography_scrna/genome_files/gene_id_feature_type.csv"),
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