# Description -------------------------------------------------------------
# Here we will prepare files to create anndata objects which can be processed in R

# Set up ------------------------------------------------------------------
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
  # Remove barcode column where there is only NA values - otherwise it will not join the data properly
  purrr::modify_depth(1, ~ .x |> dplyr::select_if(~ !all(is.na(.))))

# Divide Kang into cell and nuclei
meta_harmonized_list[["Kang_cell"]] <- meta_harmonized_list[["Kang"]] |>dplyr::filter(cell_nuclei == "cell") |> 
  dplyr::mutate(name = paste0(name, "_cell"))
meta_harmonized_list[["Kang_nuclei"]] <- meta_harmonized_list[["Kang"]] |>dplyr::filter(cell_nuclei == "nuclei") |> 
  dplyr::mutate(name = paste0(name, "_nuclei"))
meta_harmonized_list[["Kang"]] <- NULL

## Combine harmonized meta data and quality control
ann_meta_data <- purrr::imap(quality_control, function(df_qc, name){
  output <- dplyr::left_join(x = df_qc, y = meta_harmonized_list[[name]])
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

rm(quality_control, meta_harmonized)
gc()
# Add prefix to barcodes --------------------------------------------------

# Start redirecting console output to a log file
sink("my_log_file.txt")

purrr::iwalk(ann_meta_data, function(meta, name) { # Iterate over all elements of ann_meta_data
  
  # Extract study-level metadata
  platform <- unique(meta$platform)
  cell_nuclei <- unique(meta$cell_nuclei)
  
  # Basic validation: ensure unique values
  if (length(platform) > 1) {
    stop("Multiple platforms found for study '", name, "'. Each study should have one platform.")
  }
  if (length(cell_nuclei) > 1) {
    stop("Multiple cell_nuclei types found for study '", name, "'. Each study should have one type.")
  }
  
  print(paste0("Processing Study: ", name))
  print(paste0("Platform: ", platform))
  print(paste0("Cell/Nuclei: ", cell_nuclei))
  
  # Specific for Kang study naming convention for file paths
  current_study_name_for_path <- name
  if (name %in% c("Kang_nuclei", "Kang_cell")) {
    current_study_name_for_path <- "Kang"
    print(paste0("Adjusted study name for path: ", current_study_name_for_path))
  }
  
  # Define base paths to search for barcodes.tsv
  search_base_dirs <- purrr::map_chr(unique(meta$sample), function(s_name) {
    base::paste0("/work/scRNAseq/", current_study_name_for_path, "/Preprocessed/", s_name, "/Solo.out/",
                 ifelse(cell_nuclei == "cell", "Gene", "GeneFull_Ex50pAS"), "/raw")
  })
  
  # Find barcode.tsv files within these specific directories
  paths <- purrr::map(search_base_dirs, ~base::list.files(.x, pattern = "barcodes.tsv", recursive = FALSE, full.names = TRUE)) %>%
    base::unlist()
  
  # Check if any paths were found.
  if (length(paths) == 0) {
    warning("No barcode files found for study: ", name, ". Skipping this study.")
    return(NULL)
  }
  
  # Name the paths with the original sample names extracted robustly
  original_sample_names_from_paths <- purrr::map_chr(paths, ~ {
    extracted_name <- stringr::str_extract(.x, "(?<=Preprocessed/)[^/]+")
    if (is.na(extracted_name)) {
      warning("Failed to extract sample name from path: ", .x)
    }
    return(extracted_name)
  })
  
  # Add names to paths
  base::names(paths) <- original_sample_names_from_paths
  
  # Filter out any paths where name extraction failed
  paths <- paths[!is.na(base::names(paths))]
  if (length(paths) == 0) {
    warning("All barcode file paths had unextractable sample names for study: ", name, ". Skipping this study.")
    return(NULL)
  }
  
  print(paste0("Total number of barcode files found for this study: ", length(paths)))

  # Chunking Setup ----
  chunk_size <- 100
  
  # Are there less or more than 100p aths ? 
  if (length(paths) <= chunk_size) {
    # If number of paths is less than or equal to chunk_size, process all as one chunk
    print(base::paste0("Number of files (", length(paths), ") is less than or equal to chunk size (", chunk_size, "). Processing as a single chunk."))
    paths_chunks <- base::list(paths) # Wrap all paths into a single list element to mimic a chunk
    # Give it a "chunk_group" name for consistent processing (e.g., "1")
    base::names(paths_chunks) <- "1"
  } else {
    # Original chunking logic for more than chunk_size files
    # Create a grouping factor for chunking
    chunk_groups <- base::cut(base::seq_along(paths), breaks = seq(0, length(paths), by = chunk_size),
                        labels = FALSE, include.lowest = TRUE)
    # Split the named 'paths' vector into a list of named 'paths' chunks
    paths_chunks <- base::split(paths, chunk_groups)
  }

  print(paste0("Processing files in ", length(paths_chunks), " chunks of ", chunk_size, " files."))
  
  # Core Processing Function for a Single Chunk ----
  process_barcode_chunk <- function(current_chunk_paths, chunk_idx) {
    print(paste0("  Processing chunk ", chunk_idx, "/", base::length(paths_chunks), "... (", base::length(current_chunk_paths), " files)"))
    
    # Read barcode.tsv files for the current chunk
    barcodes_chunk <- purrr::map(current_chunk_paths, vroom::vroom_lines)
    
    # Check for any names in 'barcodes_chunk' that are not in 'meta$sample'
    if (!base::all(base::names(barcodes_chunk) %in% meta$sample)) {
      mismatched_names <- base::names(barcodes_chunk)[!base::names(barcodes_chunk) %in% meta$sample]
      warning("Mismatch: The following extracted path names in chunk ", chunk_idx, " are not found in meta$sample for study '", name, "': ",
              paste(mismatched_names, collapse = ", "))
      # Filter them out.
      barcodes_chunk <- barcodes_chunk[base::names(barcodes_chunk) %in% meta$sample]
      if (length(barcodes_chunk) == 0) {
        warning("No valid barcodes left in chunk ", chunk_idx, " after meta data mapping check for study: ", name, ". Skipping this chunk.")
        return(NULL) # Return NULL to indicate this chunk was skipped
      }
    }
    
    # Create mapping vectors (defined outside to be efficient, but used inside)
    sample_to_ic_id_map <- setNames(meta$ic_id_sample, meta$sample)
    ic_id_to_sample_map <- setNames(meta$sample, meta$ic_id_sample)
    
    # Rename barcode list names from original sample names to ic_id_sample for this chunk
    # So we'll first rename the list elements
    base::names(barcodes_chunk) <- sample_to_ic_id_map[base::names(barcodes_chunk)]
    
    if (any(is.na(base::names(barcodes_chunk)))) {
      warning("Some barcode names failed to map to ic_id_sample names during first rename in chunk ", chunk_idx, " for study '", name, "'. Check meta data and original names.")
      barcodes_chunk <- barcodes_chunk[!is.na(base::names(barcodes_chunk))]
      if (length(barcodes_chunk) == 0) {
        warning("No barcodes left after filtering unmapped names in chunk ", chunk_idx, " for study: ", name, ". Skipping this chunk.")
        return(NULL)
      }
    }
    
    # Barcode content manipulation based on platform for this chunk ---
    if (platform %in% c("droplet", "plate_barcode")) {
      barcodes_chunk_processed <- purrr::imap(barcodes_chunk, function(string_vector, current_ic_id_name) {
        base::paste0(current_ic_id_name, "_", string_vector)
      })
      
      print(paste0("  Barcode content prefixed for '", platform, "' platform in chunk ", chunk_idx, ".")) 
      
    } else if (platform %in% c("plate")) {
      barcodes_chunk_processed <- purrr::imap(barcodes_chunk, function(string_vector, current_ic_id_name) {
        return(current_ic_id_name)
      })
      print(paste0("  Barcode content replaced with ic_id_sample name for '", platform, "' platform in chunk ", chunk_idx, ".")) 
    } else {
      # If platform is neither, do not modify content
      barcodes_chunk_processed <- barcodes_chunk
      print(paste0("  Barcode content not modified for platform '", platform, "' in chunk ", chunk_idx, "."))
    }
    
    # End barcode content manipulation for this chunk ---
    
    
    # Rename barcode list names back from ic_id_sample to original sample names (for saving)
    base::names(barcodes_chunk_processed) <- ic_id_to_sample_map[names(barcodes_chunk_processed)]
    
    if (any(is.na(base::names(barcodes_chunk_processed)))) {
      warning("Some barcode names failed to map back to original sample names during reverse rename in chunk ", chunk_idx, " for study '", name, "'. Check meta data and ic_id names.")
      barcodes_chunk_processed <- barcodes_chunk_processed[!is.na(base::names(barcodes_chunk_processed))]
      if (length(barcodes_chunk_processed) == 0) {
        warning("No barcodes left after filtering unmapped names for saving in chunk ", chunk_idx, " for study: ", name, ". Skipping this chunk.")
        return(NULL)
      }
    }
    
    # Save the updated barcodes for the current chunk
    purrr::iwalk(barcodes_chunk_processed, function(file_content, sample_name_for_path) {
      output_dir_base <- base::paste0("/work/scRNAseq/", current_study_name_for_path, "/Preprocessed/", sample_name_for_path, "/Solo.out/")
      
      output_final_dir <- if (cell_nuclei == "cell") {
        base::paste0(output_dir_base, "Gene/raw/")
      } else if (cell_nuclei == "nuclei") {
        base::paste0(output_dir_base, "GeneFull_Ex50pAS/raw/")
      }

      output_path <- base::paste0(output_final_dir, "barcodes_prefixed.tsv")
      
      print(paste0("    Attempting to write barcode file to: ", output_path)) 
      vroom::vroom_write_lines(file_content, output_path)
      print(paste0("    Successfully wrote to: ", output_path))
    })
    
    # No need for explicit rm() and gc() here, as purrr will manage scope
    # The chunked processing naturally helps with memory
    print(paste0("  Finished chunk ", chunk_idx, "."))
  }
  
  # Apply the processing function to each chunk using purrr::iwalk ---
  # purrr::iwalk is used here because the primary effect is the side-effect of saving files
  purrr::iwalk(paths_chunks, process_barcode_chunk)
  
})

# Stop redirecting console output
sink()

# path = 'immune/raw/'
# path_out = 'immune/anndata/'
# dataset = 'GSE129785_scATAC-Hematopoiesis-All'
# X = sc.read_mtx(f'{path}{dataset}.mtx', )
# X = X.X.T
# X = sc.AnnData(X,
#                obs = pd.read_csv(f'{path}{dataset}.metadata.txt',sep = '\t', index_col = -1),
#                var = pd.read_csv(f'{path}{dataset}.features.txt',sep = '\t', index_col = 0))
# X.obsm['X_umap'] = X.obs[['UMAP1','UMAP2']].values
# X.write_h5ad(f'{path_out}{dataset}.h5ad')
