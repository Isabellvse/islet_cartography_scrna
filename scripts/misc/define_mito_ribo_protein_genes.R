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

