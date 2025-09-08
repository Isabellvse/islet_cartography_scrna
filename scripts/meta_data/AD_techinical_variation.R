# Description -------------------------------------------------------------
# Extract meta data which will be used to remove technical variance during integration
# will be used to create the first h5ad objects
# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/metadata_harmonized/"))
set.seed(1000)

# Load --------------------------------------------------------------------
## list of IC IDS ----
icid_list <- qs2::qs_read(here::here("islet_cartography_scrna/data/metadata/id_list.qs2"))

## Meta data as a list ----
meta_list <- qs2::qs_read(here::here("islet_cartography_scrna/data/metadata/meta_list_1.qs2"))

## Global study variance ----
global_variance <- vroom::vroom(here::here("islet_cartography_scrna/data/global_study_variance.csv"), 
                                delim = ";",
                                col_names = TRUE)

## Cell type annotation ----
# Fang
fang_cell_type <- readxl::read_xlsx(here::here("islet_cartography_scrna/data_raw/meta/Fang/Fang_cell_type_annotation.xlsx"),
                                     sheet = "Celltype.info", skip = 2) |> 
  dplyr::mutate(donor = paste0("Dropseq-", Donor),
                sample = donor)

# Shrestha
shrestha_cell_type <- readr::read_delim(here::here("islet_cartography_scrna/data_raw/meta/Shrestha/Shrestha.txt"), 
                                        skip = 1, 
                                        col_names = c("barcode", "Library_ID", "Age", "Sex", "Batch", "Replicates", "CellTypes")) |> 
  dplyr::mutate(Age = stringr::str_remove(Age, "y") |> as.numeric()) |> 
  dplyr::filter(Age >= 20) |> 
  dplyr::select(barcode, Library_ID, CellTypes)

# Motakis 
# motakis_seurat <- anndata::read_h5ad(here::here("islet_cartography_scrna/Motakis.h5ad"))
# motakis_seurat$obs|>
#   tibble::rownames_to_column("barcode") |>
#   vroom::vroom_write(
#     here::here(
#       "islet_cartography_scrna/data_raw/meta/Motakis/Motakis_meta_seurat.csv"
#     ),
#     delim = ";"
#   )

motakis_cell_type <- vroom::vroom(here::here("islet_cartography_scrna/data_raw/meta/Motakis/Motakis_meta_seurat.csv"), 
                               delim = ";",
                               col_names = TRUE) |> 
  dplyr::select(barcode, inferred_cell_type = cell_type, donor = Islet) |> 
  dplyr::mutate(inferred_cell_type = dplyr::case_when(
    inferred_cell_type == "endothelial cell" ~ "endothelial",
    inferred_cell_type == "leukocyte" ~ "leukocyte",
    inferred_cell_type == "mast cell" ~ "mast",
    inferred_cell_type == "pancreatic A cell" ~ "alpha",
    inferred_cell_type == "pancreatic acinar cell" ~ "acinar",
    inferred_cell_type == "pancreatic D cell" ~ "delta",
    inferred_cell_type == "pancreatic ductal cell" ~ "ductal",
    inferred_cell_type == "pancreatic epsilon cell" ~ "epsilon",
    inferred_cell_type == "pancreatic PP cell" ~ "gamma",
    inferred_cell_type == "pancreatic stellate cell" ~ "stellate",
    inferred_cell_type == "Schwann cell" ~ "Schwann",
    inferred_cell_type == "type B pancreatic cell" ~ "beta"))

# Tritschler
# tritchler <- anndata::read_h5ad(here::here("islet_cartography_scrna/Tritschler.h5ad"))
# tritchler$obs |> 
#     tibble::rownames_to_column("barcode") |>
#     vroom::vroom_write(
#       here::here(
#         "islet_cartography_scrna/data_raw/meta/Tritschler/Tritschler_meta_h5ad.csv"
#       ),
#       delim = ";"
#     )
tritchler_cell_type <- vroom::vroom(here::here("islet_cartography_scrna/data_raw/meta/Tritschler/Tritschler_meta_h5ad.csv"), 
                                  delim = ";",
                                  col_names = TRUE) |> 
  dplyr::select(barcode, donor = id, inferred_cell_type = louvain_anno_broad)


## Missing meta data ----
# Mauvais
# Loaded the seurat object from GEO and saved as a csv file - removed again as it is large
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE266291&format=file&file=GSE266291%5Fprocessed%5Frna%2Eqs%2Egz
# mauvais_Jarvis <- qs::qread(here::here("islet_cartography_scrna/Mauvais_Jarvis.qs"))
# mauvais_Jarvis@meta.data |>
#   tibble::rownames_to_column("barcode") |>
#   dplyr::filter(!grepl("HPAP-", Library)) |>
#   vroom::vroom_write(
#     here::here(
#       "islet_cartography_scrna/data_raw/meta/Mauvais_Jarvis/Mauvais_Jarvis_meta_seurat.csv"
#     ),
#     delim = ";"
#   )

mj_seurat_meta <- vroom::vroom(here::here("islet_cartography_scrna/data_raw/meta/Mauvais_Jarvis/Mauvais_Jarvis_meta_seurat.csv"), 
                                         delim = ";",
                                         col_names = TRUE) |> 
  dplyr::select(Library, age_years = age, bmi, ethnicity = ancestry) |> 
  dplyr::distinct()

## Raw meta data ----
mj_meta <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Mauvais_Jarvis/GSE266291_series_matrix.txt"), getGPL = F)
mj_meta <- Biobase::phenoData(mj_meta)@data |> 
  dplyr::mutate(donor = title,
                donor = stringr::str_remove(donor, "Pancreatic islets, non-diabetic, biological replicate "),
                donor = base::paste0("Islet_r", donor))

# Preprocess raw metadata -------------------------------------------------
# Data such as age, sex, bmi, and ethnicity was not stated in the series matrix from Mauvais Jarvis, but it is in the Seurat object

# Get ids that match between the ones I use, and the ones that is used in the seurat object meta data
mj_meta <- mj_meta |> dplyr::mutate(Library = stringr::str_extract_all(supplementary_file_1, 
                                                            stringr::str_flatten(unique(mj_seurat_meta$Library), collapse = "|")) |> 
                           as.character()) |> 
  dplyr::select(donor, Library) |> 
  dplyr::left_join(y = mj_seurat_meta) |> 
  dplyr::select(-Library)

# Preprocess cell type annotation -----------------------------------------
# Preprocess cell type annotations from studies
fang_study_annotation <- icid_list[["Fang"]] |> 
  dplyr::left_join(y = fang_cell_type, relationship = "one-to-many") |> 
  dplyr::mutate(barcode = paste0(ic_id, "_", stringr::str_remove(Cell.ID, "_(.*)")),
                inferred_cell_type = dplyr::case_when(celltype == "NA" ~ NA,
                                                    .default =  snakecase::to_snake_case(celltype))) |> 
  dplyr::select(name, sample, donor, barcode, inferred_cell_type)

shrestha_study_annotation <- meta_list[["Shrestha"]] |> 
  dplyr::select(name, donor, sample, Library_ID = sample_id) |> 
  dplyr::left_join(y = shrestha_cell_type, relationship = "one-to-many") |> 
  dplyr::left_join(y = icid_list[["Shrestha"]]) |> 
  dplyr::mutate(barcode = paste0(ic_id, "_", stringr::str_remove(barcode, "-(.*)")),
                inferred_cell_type = snakecase::to_snake_case(CellTypes)) |> 
  dplyr::select(name, sample, donor, barcode, inferred_cell_type)

# In motakis, barcodes can have multiple celltype annotations, I will combine these annotations 
# into one, so that I do not have duplicated values
motakis_study_annotation <- icid_list[["Motakis"]] |> 
  dplyr::left_join(y = motakis_cell_type, relationship = "one-to-many") |> 
  dplyr::mutate(barcode = paste0(ic_id, "_", stringr::str_remove(barcode, ".*_"))) |> 
  dplyr::select(name, sample, donor, barcode, inferred_cell_type) |> 
  dplyr::group_by(dplyr::pick(-inferred_cell_type)) |> 
  dplyr::summarise(inferred_cell_type = base::paste0(BiocGenerics::unique(inferred_cell_type), 
                                                     collapse = "_"),
                   .groups = "drop") 

tritchler_study_annotation <- icid_list[["Tritschler"]] |> 
  dplyr::left_join(y = tritchler_cell_type, relationship = "one-to-many") |> 
  dplyr::mutate(barcode = paste0(ic_id, "_", stringr::str_remove(barcode, "-(.*)"))) |> 
  dplyr::select(name, donor, barcode, inferred_cell_type)

# Preprocess --------------------------------------------------------------
## Select columns to keep ----
meta_list_sub <- purrr::map(meta_list, ~dplyr::select(.x, tidyselect::any_of(meta_cols))) 

# We use the mean time when ranges are given
## Adjust columns
meta_list_sub_cols <- meta_list_sub |> 
  purrr::map(~dplyr::rename_with(.x, .fn = snakecase::to_snake_case)) |> 
  purrr::map_at(.at = "Avrahami", ~ dplyr::mutate(.x, 
                                                  islet_center = "integrated_islet_distribution_program_iidp",
                                                  tissue = "islet",
                                                  islet_culture_medium = "prodo_islet_media",
                                                  islet_culture_medium_glucose_milimolar = 5.8
                                                  )) |> 
  purrr::map_at(.at = "Baron", ~ dplyr::mutate(.x,
                                               islet_center = "the_national_disease_research_interchange_ndri",
                                               tissue = "islet",
                                               islet_culture_medium = "cmrls",
                                               islet_culture_hours = (24+48) / 2
  )) |> 
  # In camunas, donor R257 is denoted as both T2D and elevated hba1c % in the raw data
  # (camunas_cry and camunas_facs)
  # camunas_facs %>% dplyr::filter(donor == "R257") %>% dplyr::pull("disease") %>% unique()
  # camunas_cry %>% dplyr::filter(donor == "R257") %>% dplyr::pull("disease") %>% unique()
  # I will combine these annotations into one
  purrr::map_at(.at = "Camunas", ~ dplyr::mutate(.x, 
                                                 plate = as.character(plate),
                                                 islet_culture_hours = culture_time_day*24,
                                                 islet_center = "albeta_diabetes_institute_islet_core",
                                                 islet_fresh_frozen = dplyr::case_when(cryopreserved == "Yes" ~ "frozen",
                                                                                       cryopreserved == "No" ~ "fresh"),
                                                 treatment_patch = dplyr::case_when(patched == "Yes" ~ "patch",
                                                                                    patched == "No" ~ "no_patch",
                                                                                    patched == "FACS" ~ "no_patch"),
                                                 treatment_facs= dplyr::case_when(facs_sorted == "Yes" ~ "facs",
                                                                                    facs_sorted == "No" ~ "no_facs"),
                                                 dissociation_tool = dplyr::case_when(donor %in% c(paste0("R", seq(229,242)), paste0("cryo_R", seq(229,242))) ~ "hanks_based_cell_dissociation_buffer",
                                                                                      donor %in% c(paste0("R", seq(243,269) ), paste0("cryo_R", seq(243,269))) ~ "stempro_accutase"),
                                                 dissociation_method = dplyr::case_when(donor %in% c(paste0("R", seq(229,242)), paste0("cryo_R", seq(229,242))) ~ "enzyme_free",
                                                                                          donor %in% c(paste0("R", seq(243,269)), paste0("cryo_R", seq(243,269))) ~ "enzymatic"),
                                                 tissue = "islet",
                                                 islet_culture_medium = "dmem",
                                                 islet_culture_medium_glucose_milimolar = 5.5) |> 
                  dplyr::rename(cold_ischemia_hours = cold_ischemia_time_h) |> 
                  dplyr::select(-culture_time_day,
                                -facs_sorted,
                                -patched,
                                -cryopreserved,
                                -timefrom_dispersion_days) |> 
                  dplyr::rename(instrument_facs = sorting_instrument,
                                inferred_cell_type = cell_type) |> 
                  dplyr::group_by(donor)  |> 
                  dplyr::mutate(
                    disease = paste(base::unique(disease), collapse = "_")) %>%
                  dplyr::ungroup())|> 
  purrr::map_at(.at = "Dai", ~ dplyr::mutate(.x, 
                                             islet_center = snakecase::to_snake_case(islet_center),
                                             islet_fresh_frozen = snakecase::to_snake_case(fresh_or_cryo),
                                             treatment_patch = dplyr::case_when(patched_ch_1 == "Yes" ~ "patch",
                                                                                patched_ch_1 == "No" ~ "no_patch"),
                                             treatment_facs = dplyr::case_when(facs == "NO" ~ "no_facs",
                                                                               facs == "YES(non alpha/beta)" ~ "facs"),
                                             tissue = "islet",
                                             islet_culture_medium = "dmem",
                                             islet_culture_medium_glucose_milimolar = 5.5) |> 
                  dplyr::rename(inferred_cell_type = cell_type) |> 
                  dplyr::select(-fresh_or_cryo,
                                -timefrom_dispersion_days,
                                -patched_ch_1,
                                -facs) |> 
                  dplyr::rename(cold_ischemia_hours = cold_ischemia_time)) |>
  purrr::map_at(.at = "Enge", ~ dplyr::mutate(.x, 
                                              # we don't know exactly which center the islets are from 
                                              islet_center = "integrated_islet_distribution_program_iidp | ucsf_islet_isolation_core_san_francisco_ca_usa | international_institute_for_the_advancement_of_medicine_iiam" ,
                                              tissue = "islet"
  )) |> 
  purrr::map_at(.at = "Fang", ~ dplyr::mutate(.x, 
                                              islet_center = "prodo_laboratories",
                                              tissue = "islet",
                                              islet_culture_medium = "prodo_islet_media",
                                              islet_culture_medium_glucose_milimolar = 5.8,
                                              ethnicity = dplyr::case_when(
                                                donor == "Dropseq-H1" ~ "asian_filipino",
                                                donor == "Dropseq-H2" ~ "caucasian",
                                                donor == "Dropseq-H3" ~ "caucasian",
                                                donor == "Dropseq-H4" ~ "caucasian",
                                                donor == "Dropseq-H5" ~ "caucasian",
                                                donor == "Dropseq-H6" ~ "caucasian",
                                                donor == "Dropseq-T2D1" ~ "caucasian",
                                                donor == "Dropseq-T2D2" ~ "caucasian",
                                                donor == "Dropseq-T2D3" ~ "hispanic"
                                              ),
                                              hba_1_c_percent = dplyr::case_when(
                                                donor == "Dropseq-H1" ~ 5.4,
                                                donor == "Dropseq-H2" ~ 5.2,
                                                donor == "Dropseq-H3" ~ 5.0,
                                                donor == "Dropseq-H4" ~ 5.6,
                                                donor == "Dropseq-H5" ~ 4.9,
                                                donor == "Dropseq-H6" ~ 5.4,
                                                donor == "Dropseq-T2D1" ~ 8.9,
                                                donor == "Dropseq-T2D2" ~ 5.2,
                                                donor == "Dropseq-T2D3" ~ 7.1
                                              ),
                                              cause_of_death = dplyr::case_when(
                                                donor == "Dropseq-H1" ~ "stroke",
                                                donor == "Dropseq-H2" ~ "trauma",
                                                donor == "Dropseq-H3" ~ "anoxia",
                                                donor == "Dropseq-H4" ~ "stroke",
                                                donor == "Dropseq-H5" ~ "stroke",
                                                donor == "Dropseq-H6" ~ "automobile_accident",
                                                donor == "Dropseq-T2D1" ~ "anoxia",
                                                donor == "Dropseq-T2D2" ~ "cerebral_vascular_accident",
                                                donor == "Dropseq-T2D3" ~ "stroke")) |> 
                  dplyr::right_join(y = fang_study_annotation)) |> 
  purrr::map_at(.at = "Gurp", ~ dplyr::mutate(.x, 
                                              islet_center = "integrated_islet_distribution_program_iidp",
                                              instrument_facs = "moflo_astrios",
                                              treatment_facs = "facs",
                                              tissue = "islet")) %>% 
  purrr::map_at(.at = c("HPAP_10x", "HPAP_fluidigm", "HPAP_patch_22", "HPAP_patch_23"), ~ dplyr::mutate(.x, 
                                              islet_center = "human_pancreas_analysis_program_hpap") |> 
                  dplyr::rowwise() |>
                  dplyr::mutate(cold_ischemia_hours = base::round(parse_to_hours(cold_ischemia_time_ddhhmmss), 2),
                                islet_allocation_facility = snakecase::to_snake_case(allocation_via),
                                tissue = "islet") |> 
                  dplyr::ungroup() |> 
                  dplyr::select(-cold_ischemia_time_ddhhmmss, -allocation_via)) |> 
  purrr::map_at(.at = c("HPAP_patch_22", "HPAP_patch_23"), ~ 
                  dplyr::mutate(.x, 
                                treatment_patch = "patch",
                                tissue = "islet")) |> 
  purrr::map_at(.at = "HPAP_fluidigm", ~ dplyr::mutate(.x, 
                                                       islet_culture_medium = "dmem",
                                                       islet_culture_medium_glucose_milimolar = 5.5,
                                                       islet_culture_hours = (24+(24*4))/2)) |> 
  purrr::map_at(.at = "Kang", ~ 
                  dplyr::mutate(.x, 
                                islet_center = "prodo_laboratories",
                                tissue = "islet",
                                islet_culture_medium = "rpmi",
                                islet_culture_medium_glucose_milimolar = 5) |> 
                  dplyr::select(-tissue_source) |> 
                  dplyr::rename(islet_culture_hours = culture_time_hours)) |> 
  purrr::map_at(.at = "Lawlor", ~ 
                  dplyr::mutate(.x, 
                                islet_center = dplyr::case_when(donor == "ACEK420A" ~ "integrated_islet_distribution_program_iidp",
                                                                .default = "prodo_laboratories"),
                                islet_fresh_frozen = "fresh",
                                tissue = "islet",
                                islet_culture_medium = "prodo_islet_media",
                                islet_culture_medium_glucose_milimolar = 5.8,
                                islet_culture_hours = 24)) |> 
  purrr::map_at(.at = "Mauvais_Jarvis", ~ 
                  dplyr::mutate(.x, 
                                tissue = "islet",
                                islet_culture_medium = "rpmi",
                                islet_culture_medium_glucose_milimolar = 11,
                                islet_culture_hours = 24,
                                islet_center = "prodo_laboratories | integrated_islet_distribution_program_iidp") |> # we don't know exactly which center the islets are from
                  dplyr::left_join(y = mj_meta)) |> 
  purrr::map_at(.at = "Motakis", ~ dplyr::mutate(.x, islet_allocation_facility = snakecase::to_snake_case(center_ch_1),
                                                 islet_center = snakecase::to_snake_case(islet_center),
                                                 tissue = "islet",
                                                 islet_culture_medium = "cmrl",
                                                 islet_culture_hours = 14*24) |> 
                  dplyr::select(-center_ch_1) |> 
                  dplyr::left_join(y = motakis_study_annotation)) |> 
  purrr::map_at(.at = "Muraro", ~ dplyr::mutate(.x, 
                                                treatment_facs = "facs",
                                                instrument_facs = "facs_aria_ii | facs_jazz",
                                                islet_center = snakecase::to_snake_case(islet_center),
                                                tissue = "islet",
                                                islet_culture_medium = "cmrl_066",
                                                islet_culture_medium_glucose_milimolar = 5.5,
                                                islet_culture_hours = ((3*24)+(5*24))/2)) |> 
  purrr::map_at(.at = "Segerstolpe", ~ dplyr::mutate(.x, treatment_facs = "facs",
                                                     islet_center = "prodo_laboratories",
                                                     tissue = "islet",
                                                     islet_culture_medium = "prodo_islet_media",
                                                     islet_culture_medium_glucose_milimolar = 5.8,
                                                     islet_culture_hours = 96)) %>% 
  purrr::map_at(.at = "Shrestha", ~ dplyr::rename(.x, islet_culture_hours = culture_time_hours,
                                                  cold_ischemia_hours = cold_ischaemia_time_hours) |> 
                  dplyr::mutate(islet_center = snakecase::to_snake_case(islet_center),
                                tissue = "islet",
                                islet_culture_medium = "cmrl_066",
                                islet_culture_medium_glucose_milimolar = 5.5) |> 
                  dplyr::right_join(y = shrestha_study_annotation)) |> 
  purrr::map_at(.at = "Son", ~ dplyr::mutate(.x, islet_center = "integrated_islet_distribution_program_iidp",
                                             tissue = "islet",
                                             islet_culture_medium = "prodo_islet_media",
                                             islet_culture_medium_glucose_milimolar = 5.8)) |> 
  purrr::map_at(.at = "Tritschler", ~ dplyr::mutate(.x, islet_center = "albeta_diabetes_institute_islet_core",
                                             tissue = "islet",
                                             islet_culture_medium = "cmrl_066",
                                             islet_culture_medium_glucose_milimolar = 5.5) |> 
                  dplyr::right_join(y = tritchler_study_annotation)) |> 
  purrr::map_at(.at = "Wang", ~ dplyr::mutate(.x, islet_center = snakecase::to_snake_case(islet_center),
                                                    tissue = "islet",
                                              disease = dplyr::case_when(disease == "control" ~ "nd",
                                                                         .default = as.character(base::tolower(disease))),
                                              islet_culture_hours = as.numeric(cultured_days)*24,
                                              islet_culture_medium = "prodo_islet_media",
                                              islet_culture_medium_glucose_milimolar = 5.8,
                                              disease = stringr::str_to_lower(disease)) |> 
                  dplyr::select(-cultured_days)) |> 
  purrr::map_at(.at = "Wang_Sander", ~ dplyr::mutate(.x, islet_center = snakecase::to_snake_case(islet_center),
                                                    tissue = "islet",
                                                    islet_fresh_frozen = "frozen")) |> 
  purrr::map_at(.at = "Xin", ~ dplyr::mutate(.x, islet_center = "prodo_laboratories",
                                                     tissue = "islet")) |> 
  purrr::map_at(.at = "Xin_Diabetes", ~ dplyr::mutate(.x, islet_center = "prodo_laboratories",
                                             tissue = "islet")) |> 
  purrr::map_at(.at = "Zhang", ~ dplyr::mutate(.x, islet_center = "prodo_laboratories",
                                                      tissue = "islet",
                                               islet_culture_medium = "prodo_islet_media",
                                               islet_culture_medium_glucose_milimolar = 5.8,
                                               islet_culture_hours = 16)) # they specify overnight

# Combining with additional data ------------------------------------------
## Global variance + some more processing ----
meta_sub_cols <- meta_list_sub_cols |> 
  purrr::list_rbind() |> 
  dplyr::mutate(instrument_seq = snakecase::to_snake_case(instrument_model)) |>
  dplyr::left_join(y = global_variance, relationship = "many-to-one", by = "name") |> # Combine with global variance
  # Only add values from y if x is NA
  dplyr::mutate(dissociation_tool = dplyr::coalesce(dissociation_tool.x, dissociation_tool.y),
                # if we don't know exactly, we give it NA
                dissociation_tool = dplyr::case_when(dissociation_tool == "enzymatic|Accutase|enzymefree_buffer|Hanks’-based Cell Dissociation Buffer" ~ NA, 
                                                     .default = as.character(dissociation_tool)),
                # Only add values from y if x is NA
                dissociation_method = dplyr::coalesce(dissociation_method.x, dissociation_method.y),
                # if we don't know exactly, we give it NA
                dissociation_method = dplyr::case_when(dissociation_method == "enzymatic|nonenzymatic" ~ NA,
                                                       .default = as.character(dissociation_method)),
                # Only add values from y if x is NA
                instrument_seq = dplyr::coalesce(instrument_seq.x, instrument_seq.y)) |> 
  # Remove columns
  dplyr::select(-tidyselect::ends_with(".y"), 
                -tidyselect::ends_with(".x"), 
                -instrument_model,
                - source_tissue, 
                study_cell_annotation = inferred_cell_type) # rename inferred_cell_type to study_cell_annotation

## Islet cartography ids ----
meta_sub_id <- icid_list |> 
  purrr::list_rbind() |> 
  dplyr::select(-ic_id) |> # remove ic_id column
  dplyr::mutate(ic_id_donor = paste0(ic_id_study, "_", ic_id_donor), # combine study and donor id
                ic_id_sample = paste0(ic_id_donor, "_", ic_id_sample), # combine study, donor, and sample id
                dplyr::across(.cols = tidyselect::starts_with("ic_id"), .fns = ~paste0("ic", "_", .x))) |> #add ic_ prefix
  dplyr::left_join(y = meta_sub_cols) # Join with meta data

## Order of columns - check this 
meta_harmonized  <- meta_sub_id |> dplyr::select(tidyselect::all_of(meta_variance_ordered))

# Harmonize columns -------------------------------------------------------

## Study annotation ----
meta_harmonized <- meta_harmonized %>%
  dplyr::mutate(
    study_cell_annotation_lower = stringr::str_to_lower(study_cell_annotation),
    study_cell_annotation_harmonized = dplyr::case_when(
      
      # NA handling
      base::is.na(study_cell_annotation) ~ "unknown",
      
      # Beta cells
      study_cell_annotation_lower %in% c("beta", "beta cell", "all_beta") ~ "beta",
      
      # Alpha cells
      study_cell_annotation_lower %in% c("alpha", "alpha cell") ~ "alpha",
      
      # Delta cells
      study_cell_annotation_lower %in% c("delta", "delta cell") ~ "delta",
      
      # Gamma / PP cells
      study_cell_annotation_lower %in% c("gamma", "pp", "gamma/pp", "gamma cell") ~ "gamma",
      
      # Epsilon
      study_cell_annotation_lower %in% c("epsilon", "epsilon cell") ~ "epsilon",
      
      # Ductal cells
      study_cell_annotation_lower %in% c("duct", "ductal", "ductal cell") ~ "ductal",
      
      # Acinar cells
      study_cell_annotation_lower %in% c("acinar", "acinar cell") ~ "acinar",
      
      # Mesenchymal cells
      study_cell_annotation_lower %in% c("mesenchyme", "mesenchymal") ~ "mesenchymal",
      
      # Stellate / PSC cells
      study_cell_annotation_lower %in% c("stellate", "psc", "psc cell") ~ "stellate",
      
      # Endothelial cells
      study_cell_annotation_lower %in% c("endothelial", "endothelial cell") ~ "endothelial",
      
      # Immune-related
      study_cell_annotation_lower == "leukocyte" ~ "leukocyte",
      study_cell_annotation_lower %in% c("mast", "mast cell") ~ "mast",
      study_cell_annotation_lower == "immune" ~ "immune",
      
      # Schwann
      study_cell_annotation_lower == "schwann" ~ "schwann",
      
      # QC-related or uncertain
      study_cell_annotation_lower %in% c("fail_qc", "lost", "dropped", "masked", "not applicable") ~ "excluded",
      
      # Unclassified or other
      study_cell_annotation_lower %in% c(
        "unclassified cell", 
        "unclassified endocrine cell", 
        "unsure", 
        "none/other", 
        "other"
      ) ~ "unknown",
      
      study_cell_annotation_lower == "mhc class ii cell" ~ "mhc_class_ii",
      study_cell_annotation_lower == "co-expression cell" ~ "co_expression",
      
      # Any label containing two cell types joined with "_"
      stringr::str_detect(study_cell_annotation_lower, ".*_.*") ~ "co_expression"
    )
  ) %>%
  dplyr::select(-study_cell_annotation_lower)

## Ethnicity borad categories ----
meta_harmonized  <- meta_harmonized  %>%
  dplyr::mutate(ethnicity_broad_harmonized = dplyr::case_when(
    # Of European Descent
    stringr::str_to_lower(ethnicity) %in% c("european_american", "caucasian", "white") ~ "of_european_descent",
    
    # Of African Descent
    stringr::str_to_lower(ethnicity) %in% c("african_american", "black", "africanamerican", "black_or_african_american") ~ "of_african_descent",
    
    # Of Latin American / Hispanic Descent
    stringr::str_to_lower(ethnicity) %in% c("hispanic", "hispanic_latino") ~ "of_latin_american_hispanic_descent",
    
    # Of Asian Descent (Broad)
    stringr::str_to_lower(ethnicity) %in% c("asian", "asian_filipino", "asian_indian") ~ "of_asian_descent",
    
    # Unknown / Not Specified (handles both "na" string and actual NA)
    is.na(ethnicity) | stringr::str_to_lower(ethnicity) == "na" ~ NA,
  ))

## Ethnicity sub categories ----
meta_harmonized  <- meta_harmonized  %>%
  dplyr::mutate(ethnicity_sub_harmonized = dplyr::case_when(
    # Unknown / Not Specified (handles both "na" string and actual NA)
    is.na(ethnicity) | stringr::str_to_lower(ethnicity) == "na" ~ NA,
    
    # Of European Descent
    stringr::str_to_lower(ethnicity) %in% c("european_american", "caucasian", "white") ~ "of_european_descent",
    
    # Of African Descent
    stringr::str_to_lower(ethnicity) %in% c("african_american", "black", "africanamerican", "black_or_african_american") ~ "of_african_descent",
    
    # Of Latin American Descent
    stringr::str_to_lower(ethnicity) == "hispanic_latino" ~ "of_latin_american_descent",
    stringr::str_to_lower(ethnicity) == "hispanic" ~ "of_hispanic_descent",
    
    # Of Asian Descent (Broad)
    stringr::str_to_lower(ethnicity) == "asian" ~ "of_asian_descent",
    stringr::str_to_lower(ethnicity) == "asian_filipino" ~ "of_filipino_descent",
    stringr::str_to_lower(ethnicity) == "asian_indian" ~ "of_indian_descent",
  ))

## Cause of death ----
meta_harmonized <- meta_harmonized %>%
  dplyr::mutate(
    cause_of_death_broad_harmonized = dplyr::case_when(

      # Anoxia (covers general anoxia, including asthma-related)
      stringr::str_detect(stringr::str_to_lower(cause_of_death), "anoxia") ~ "anoxia",
      
      # Cerebrovascular Event / Stroke (includes stroke, CVA, aneurysm)
      stringr::str_detect(stringr::str_to_lower(cause_of_death), "stroke") |
        stringr::str_detect(stringr::str_to_lower(cause_of_death), "cerebrovascular") |
        stringr::str_to_lower(cause_of_death) == "aneurysm" ~ "cerebrovascular_disease",
      cause_of_death == "cerebral_vascular_accident" ~ "cerebrovascular_disease",
      
      # Head Trauma (covers specific head traumas and related incidents like MVA or fall)
      stringr::str_detect(stringr::str_to_lower(cause_of_death), "head_trauma|head trauma") ~ "head_trauma",
      
      # Accident (Motor Vehicle) - for standalone auto accidents not specified as head trauma
      stringr::str_to_lower(cause_of_death) == "automobile_accident" ~ "accident_motor_vehicle",
      
      # Trauma (Other / Unspecified) - for general "trauma"
      stringr::str_to_lower(cause_of_death) == "trauma" ~ "trauma_unspecified",
      
      # NR - not registered
      cause_of_death == "NR" ~ NA_character_),
    
    cause_of_death_sub_harmonized = dplyr::case_when(
      # Drug Intoxication (specific pattern, check before more general 'anoxia')
      stringr::str_detect(stringr::str_to_lower(cause_of_death), "drug_intoxication") ~ "drug_intoxication",
      
      stringr::str_detect(stringr::str_to_lower(cause_of_death), "stroke|cerebral_vascular_accident") ~ "stroke",
      cause_of_death == "Head trauma - gunshot wound/suicide" ~ "gunshot",
      cause_of_death == "anoxia_asthma" ~ "asthma",
      cause_of_death == "head_trauma|fall" ~ "fall",
      .default = as.character(cause_of_death_broad_harmonized)
    )
  )
# Save --------------------------------------------------------------------
# as qs object
qs2::qs_save(meta_harmonized, here::here("islet_cartography_scrna/data/metadata_harmonized/meta_harmonized.qs2"))
# as csv file
vroom::vroom_write(meta_harmonized, here::here("islet_cartography_scrna/data/metadata_harmonized/meta_harmonized.csv"), delim = ",", col_names = TRUE)


# meta_harmonized_list <- meta_harmonized |>
#   # Split by study name
#   (\(df) base::split(df, factor(df$name)))() |>
#   # Remove columns where there is only NA values - otherwise it will not join the data properly
#   purrr::modify_depth(1, ~ .x |> dplyr::select_if( ~ !all(is.na(.))))
# 
# # Divide Kang into cell and nuclei
# meta_harmonized_list[["Kang_cell"]] <- meta_harmonized_list[["Kang"]] |>
#   dplyr::filter(cell_nuclei == "cell") |>
#   dplyr::mutate(name = paste0(name, "_cell"))
# meta_harmonized_list[["Kang_nuclei"]] <- meta_harmonized_list[["Kang"]] |>
#   dplyr::filter(cell_nuclei == "nuclei") |>
#   dplyr::mutate(name = paste0(name, "_nuclei"))
# meta_harmonized_list[["Kang"]] <- NULL
# 
# meta_harmonized_old <- qs2::qs_read(
#   here::here(
#     "islet_cartography_scrna/data/metadata_harmonized/meta_harmonized.qs2"
#   )
# )
# 
# 
# meta_harmonized_old_list <- meta_harmonized_old |>
#   # Split by study name
#   (\(df) base::split(df, factor(df$name)))() |>
#   # Remove columns where there is only NA values - otherwise it will not join the data properly
#   purrr::modify_depth(1, ~ .x |> dplyr::select_if( ~ !all(is.na(.))))
# 
# # Divide Kang into cell and nuclei
# meta_harmonized_old_list[["Kang_cell"]] <- meta_harmonized_old_list[["Kang"]] |>
#   dplyr::filter(cell_nuclei == "cell") |>
#   dplyr::mutate(name = paste0(name, "_cell"))
# meta_harmonized_old_list[["Kang_nuclei"]] <- meta_harmonized_old_list[["Kang"]] |>
#   dplyr::filter(cell_nuclei == "nuclei") |>
#   dplyr::mutate(name = paste0(name, "_nuclei"))
# meta_harmonized_old_list[["Kang"]] <- NULL
# 
# all.equal(names(meta_harmonized), names(meta_harmonized_old))
# map2(meta_harmonized_list, meta_harmonized_old_list, all.equal)
