# Description -------------------------------------------------------------
# Here we will prepare files to create anndata objects which can be processed in R

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(c(here::here("islet_cartography_scrna/data/anndata/raw"),
                   here::here("islet_cartography_scrna/data/anndata/metadata")))
set.seed(1000)
# Load --------------------------------------------------------------------
## Harmonized meta data ----
meta_harmonized <- qs2::qs_read(here::here("islet_cartography_scrna/data/metadata_harmonized/meta_harmonized.qs2")) 

## Quality control data ----
## Quality control ----
paths <- base::list.files(
  path = here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics_updated"),
  pattern = "quality_metrics.csv",
  full.names = TRUE,
  recursive = TRUE
)

# set names
base::names(paths) <- purrr::map_chr(paths, ~ stringr::str_extract(.x, pattern = "(?<=quality_metrics_updated/)[^/]+(?=_quality_metrics\\.csv)"))

# load
quality_control <- purrr::map(paths, ~ vroom::vroom(.x, delim = ",", col_names = TRUE) |>
                                dplyr::select(name,
                                              sample,
                                              cell_nuclei,
                                              library_prep,
                                              rna_count,
                                              barcode,
                                              # Use any_of to select nCounts and nUMIs if they exist
                                              tidyselect::any_of(c("nCounts", "nUMIs")),
                                              nFeatures,
                                              mitochondrial_fraction,
                                              coding_fraction,
                                              contrast_fraction,
                                              complexity,
                                              # Renaming directly during select for these specific cases
                                              uniquely_mapped_reads_percent = `Uniquely_mapped_reads_%`,
                                              unmapped_reads_percent = `Unmapped_reads_%`,
                                              excluded) |>
                                # Add a rename step here, using any_of to target 'nUMIs' if it's present
                                dplyr::rename(tidyselect::any_of(c(n_umi = "nUMIs",
                                                                   n_feature = "nFeatures",
                                                                   n_count = "nCounts")))
                              )

# Preprocess --------------------------------------------------------------
meta_harmonized_list <- meta_harmonized |> 
  # Split by study name
  (\(df) base::split(df, factor(df$name)))() |> 
  # Remove columns where there is only NA values - otherwise it will not join the data properly
  purrr::modify_depth(1, ~ .x |> dplyr::select_if(~ !all(is.na(.))))

# Divide Kang into cell and nuclei
meta_harmonized_list[["Kang_cell"]] <- meta_harmonized_list[["Kang"]] |>dplyr::filter(cell_nuclei == "cell") |> 
  dplyr::mutate(name = paste0(name, "_cell"))
meta_harmonized_list[["Kang_nuclei"]] <- meta_harmonized_list[["Kang"]] |>dplyr::filter(cell_nuclei == "nuclei") |> 
  dplyr::mutate(name = paste0(name, "_nuclei"))
meta_harmonized_list[["Kang"]] <- NULL

## Combine harmonized meta data and quality control
ann_meta_data <- purrr::imap(quality_control, function(df_qc, name){
  meta <- meta_harmonized_list[[name]]
  
  # Do specific subset for samples with celltype annotations, as not all barcodes in qc are in the meta data 
  if (all(c("study_cell_annotation_harmonized", "barcode") %in% colnames(meta))) {
    
    meta_sub <- meta |> 
      dplyr::select(-study_cell_annotation_harmonized, -barcode, -study_cell_annotation)|> 
      dplyr::distinct()
    
    output_tmp <- dplyr::left_join(x = df_qc, y = meta_sub)
    
    output <- output_tmp %>% dplyr::left_join(y = meta)
    
  } else {
    output <- dplyr::left_join(x = df_qc, y = meta)
  }
  # Add platform column
  output <- output |> dplyr::mutate(library_prep = base::tolower(library_prep), 
                platform = dplyr::case_when(library_prep %in% droplet_based ~ "droplet",
                                            library_prep %in% plate_based ~ "plate",
                                            library_prep %in% plate_based_bc ~ "plate_barcode"))  |>  
    dplyr::relocate(platform, .after = "library_prep")
  return(output)
})

purrr::iwalk(ann_meta_data, function(files, name){vroom::vroom_write(files, 
                                               paste0(here::here("islet_cartography_scrna/data/anndata/metadata/"), 
                                                      name, 
                                                      "_metadata.csv"), 
                                               delim = ",", 
                                               col_names = TRUE)
  })

rm(quality_control, meta_harmonized, meta_harmonized_list)
gc()

# Test a new way to add prefix --------------------------------------------
meta <- purrr::imap_dfr(ann_meta_data, ~ dplyr::select(.x, name, sample, ic_id_sample, platform, cell_nuclei)) |> 
  dplyr::mutate(name = dplyr::case_when(name %in% c("Kang_nuclei", "Kang_cell") ~ "Kang",
                                        .default = as.character(name))) %>% 
  dplyr::distinct() |> 
  dplyr::mutate(path = purrr::pmap_chr(list(name, sample, cell_nuclei), function(name, sample, cell_nuclei){
    base::paste0("/work/scRNAseq/", name, "/Preprocessed/", sample, "/Solo.out/",
                 ifelse(cell_nuclei == "cell", "Gene", "GeneFull_Ex50pAS"), "/raw/barcodes.tsv")
    
  }),
  exists = file.exists(path),
  save_path = stringr::str_replace(path, "barcodes.tsv", "barcodes_prefixed.tsv"))

sink(here::here("islet_cartography_scrna/AA_prepare_files_for_anndata.txt"))
meta |> 
  (\(df) base::split(df, factor(df$path)))() |> 
  purrr::walk(
    .f = function(df) {
      process_barcode_files(
        exists = df$exists,
        path = df$path,
        ic_id_sample = df$ic_id_sample,
        platform = df$platform,
        save_path = df$save_path
      )
    }
  )
sink()


# Check that all files have been generated --------------------------------
meta <- meta %>% dplyr::mutate(
  prefix_exists = file.exists(save_path))

table(meta$prefix_exists)
# TRUE