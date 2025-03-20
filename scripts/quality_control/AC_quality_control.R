# Description -------------------------------------------------------------
# Perform quality control on each dataset

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/quality_control/"))
set.seed(1000)
# Load --------------------------------------------------------------------
## star_quality control ----
star_quality <- qs2::qs_read(here::here("islet_cartography_scrna/data/quality_control/star_quality_raw.qs2"))

## ic ids list ----
icid_list <- qs2::qs_read(here::here("islet_cartography_scrna/data/metadata/id_list.qs2"))

## Genes for qc
genes <- qs2::qs_read(here::here("islet_cartography_scrna/data/quality_control/mito_ribo_protein_genes.qs2"))

## mtx files ----
paths <- base::paste0("/work/scRNAseq/", base::names(icid_list), "/Preprocessed/*/Solo.out")
# mtx_gene <- purrr::map(paths, ~ read_mtx(.x, feature = "Gene"))
# mtx_gene <- purrr::map2(mtx_gene, names(mtx_gene), ~ sample_to_ic_id(.x, icid_list[[.y]]))
# mtx_gene <- purrr::map(mtx_gene, prefix_colnames)
# 
# mtx_gene <- read_mtx(paths[[2]], feature = "Gene")
# mtx_gene <- sample_to_ic_id(mtx_gene, icid_list[["Baron"]])
# mtx_gene <- prefix_colnames(mtx_gene)
# 
# mtx_genefull <- read_mtx(paths[[2]], feature = "GeneFull")
# mtx_genefull <- sample_to_ic_id(mtx_genefull, icid_list[["Baron"]])
# mtx_genefull <- prefix_colnames(mtx_genefull)
# 
# barcodes <- remove_empty_dropslets(mtx_gene[["ic_2_2_2"]])
# 
# quality_met <- quality_metrics(mtx = mtx_gene[["ic_2_2_2"]],
#                 barcodes = barcodes[["not_empty"]], 
#                 mitogenes = genes[["mito_genes"]], 
#                 ribogenes = genes[["ribo_genes"]],
#                 pcgenes = genes[["protein_genes"]], 
#                 mtx_gene = mtx_gene[["ic_2_2_2"]], 
#                 mtx_genefull = mtx_genefull[["ic_2_2_2"]])

test <- process_study_samples(path = paths[[2]], study_metadata = star_quality[["Baron"]], genes = genes, study_name = "Baron", droplet_based = droplet_based, plate_based = plate_based)

test |> 
  # dplyr::select(ic_id, quality_met) %>% 
  # tibble::deframe()
# Preprocess - method -----------------------------------------------------
# Extract which sequencing method was used
star_quality |> 
  purrr::modify_depth(1, ~dplyr::select(., "name", "library_prep", "cell_nuclei") |> 
                        dplyr::distinct()) |> 
  dplyr::bind_rows() %>% 
  split(.$ic_id) 

# Obtain putative cells ---------------------------------------------------


