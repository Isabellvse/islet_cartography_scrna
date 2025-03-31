#' Create directories if they do not exist
#'
#' This function checks if each directory in the provided list exists. If any directory
#' does not exist, it creates the directory (along with any necessary parent directories).
#' It prints a message for each directory indicating whether it was created or already exists.
#'
#' @param output_dir A character vector or list of directory paths to be checked and created.
#'
#' @return NULL This function performs a side-effect (creating directories) and does not return anything.
#'
#' @export
#'
#' @examples
#' # Define a list of directories
#' output_dirs <- list("data/quality_control/rna", "data/plots", "data/other_output")
#'
#' # Call the function to create the directories
#' create_directories(output_dirs)
create_directories <- function(output_dir) {
  purrr::walk(output_dir, ~{
    if (!dir.exists(.x)) {
      dir.create(.x, recursive = TRUE)  # create the directory if it doesn't exist
      print(paste0(.x, " has been created!"))
    } else {
      print(paste0(.x, " already exists!"))
    }
  })
}


# Loading  ----------------------------------------------------------------
load_data_with_classes <- function(data_file = "data.csv", class_file = "column_classes.csv") {
  # Read the column classes
  df_classes <- vroom::vroom(class_file, delim = ",", col_types = cols(column = vroom::col_character(),
                                                                       class = vroom::col_character()))
  
  # Create a named list of column types based on the class information
  col_types <- purrr::set_names(
    purrr::map(df_classes$class, function(x) {
      base::switch(x,
                   "logical" = vroom::col_logical(),
                   "integer" = vroom::col_integer(),
                   "numeric" = vroom::col_double(),
                   "character" = vroom::col_character(),
                   "date" = vroom::col_date(format = ""),
                   "time" = vroom::col_time(format = ""),
                   "datetime" = vroom::col_datetime(format = ""),
                   "factor" = vroom::col_factor(levels = NULL, ordered = FALSE),
                   vroom::col_guess()  # Default to guessing for unknown types
      )
    }),
    df_classes$column
  )
  
  # Use vroom to load the data with the correct column types
  df_loaded <- vroom::vroom(data_file, col_types = col_types, delim = ",", n_max = 10000)
  
  return(df_loaded)
}

#' Load and Process Star Quality Data
#'
#' This function loads a data file, processes any columns containing percentage
#' signs (`%`) by removing the `%` and converting the values to numeric. It then
#' adds a `name` column, which is assigned the value of the `i` parameter, and
#' relocates this column to the first position in the resulting data frame.
#'
#' @param path A character string specifying the path to the .csv file.
#' @param i A value (numeric, character, etc.) that will be used to populate the
#'          `name` column in the data frame. This could be an identifier or a
#'          label for the data set.
#'
#' @return A data frame where:
#'   - Columns containing percentage signs (`%`) are converted to numeric by
#'     removing the `%` and converting the remaining values to numeric.
#'   - A new column named `name` is added, containing the value passed in the `i`
#'     parameter.
#'   - The `name` column is relocated to the first position.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' # df <- load_star_quality("path/to/your/file.csv", "SampleName")
#' # head(df)
#'
load_star_quality <- function(path, i){
  df <- path |> vroom::vroom()  |>
    dplyr::mutate(dplyr::across(.cols = dplyr::contains("%"),
                                .fns = ~stringr::str_remove(.x, pattern = "%") |>
                                  base::as.numeric()),
                  dplyr::across(.cols = c("Deletion_rate_per_base", "Insertion_rate_per_base"), 
                                .fns = ~stringr::str_remove(.x, pattern = "%") |>
                                  base::as.numeric()),
                  name = i) |>
    dplyr::rename("Deletion_rate_per_base_%" = "Deletion_rate_per_base",
                  "Insertion_rate_per_base_%" = "Insertion_rate_per_base") |> 
  dplyr::relocate(name)
  return(df)
}

# Quality control ---------------------------------------------------------


#' Read and Load mtx files in parallel
#'
#' @param path A string specifying the base directory containing mtx files.
#' @param feature A string specifying the folder name for feature type (e.g., "Gene" or "GeneFull")
#' @param workers An integer specifying the number of parallel workers to use. Default is the number of available CPU cores (parallelly::availableCores()).
#'
#' @returns A named list of mtx files, where names correspond to extracted sample names (from paths)
#'
#' @examples read_mtx(path = "data/Preprocessed/*/Solo.out/", feature = "Gene")
#' where * correspond to sample folder name
read_mtx <- function(path, feature){
  
  # Get list of all relevant files
  base_path <- base::paste0(path, "/", feature, "/", "raw")
  
  base::message(base::paste0("Reading mtx files from: ", base_path))
  # Get file paths
  files <- base::list.files(path = base::Sys.glob(base_path),
                            pattern = "matrix.mtx|barcodes.tsv|features.tsv",
                            recursive = TRUE,
                            full.names = TRUE)
  
  # Extract sample names from file paths
  file_info <- base::data.frame(
    path = files,
    sample = stringr::str_extract(files, "(?<=Preprocessed/)[^/]+"),
    matrix = stringr::str_extract(files, "barcodes.tsv|features.tsv|matrix.mtx"))
  
  # Reorder dataframe to ensure correct matching of paths
  file_group <- file_info |>
    tidyr::pivot_wider(names_from = matrix, values_from = path) |>
    dplyr::rename(mtx = `matrix.mtx`,
                  barcodes = `barcodes.tsv`,
                  features = `features.tsv`) |> 
    tibble::column_to_rownames("sample")
  
  mtx_list <- purrr::pmap(file_group, function(mtx, barcodes, features) {
    Seurat::ReadMtx(mtx = mtx, cells = barcodes, features = features, feature.column = 1)
  }
  )|> 
    purrr::set_names(nm = BiocGenerics::rownames(file_group)) 
  
  base::message(base::paste0("Finished reading mtx files from: ", base_path))
  
  return(mtx_list)
}

sample_to_ic_id <- function(mtx_list, ic_id_df) {
  base::names(mtx_list) <- stringi::stri_replace_all_regex(
    base::names(mtx_list),
    pattern = ic_id_df$sample,
    replacement = ic_id_df$ic_id,
    vectorize = FALSE
  )
  return(mtx_list)
}


prefix_colnames <- function(mtx_list) {
  
  purrr::imap(mtx_list, ~ {
    # Get current column names
    colnames_x <- BiocGenerics::colnames(.x)
    
    if (!base::is.null(colnames_x) && base::length(colnames_x) > 0) {
      if (base::length(colnames_x) > 1) {
        # If multiple columns exist, add a prefix
        BiocGenerics::colnames(.x) <- base::paste0(.y, "_", colnames_x)
      } else {
        # If only one column exists, replace it entirely with the sample name
        BiocGenerics::colnames(.x) <- .y
      }
    }
    
    return(.x)
  }) 
}


#' @title rank_barcodes
#'
#' @description Function to rank barcodes according to the number of UMIs (or genes) and detect a cut-off point for putative cell- (or nuclei-) containing barcodes.
#'
#' @param counts A matrix containing counts for all barcodes prior to any filtering.
#' @param type A string ("UMI" or "Genes"), which indicates which feature to use for ranking barcodes [default = "UMI"]. See details for more information.
#' @param psi.min A number indicating the lowest number of breakpoints to test when approximating the curve [default = 2]
#' @param psi.max A number indicating the highest number of breakpoints to test when approximating the curve [default = 5]
#' @param alpha A number indicating the region to find breakpoints within [default = 0.001]. See details for more information.
#' @param alpha.max The maximum allowable number for indicating the region to find breakpoints within [default = 0.05]. See details for more information.
#' @param boot A number indicating the number of bootstrap replicates used to infer breakpoints [default = 10]. Maybe necessary to increase if setting psi.max to a large number.
#' @param factor A number indicating the number of folds above the error of the best model is allowed [default = 1.5]. See details for more information.
#' @param threshold A boolean (TRUE or FALSE), which indicates if the threshold to be included in the output [default = TRUE]
#'
#' @details
#' \strong{Choosing the type of feature to use for ranking of barcodes}\cr
#' In our experience, using the number of UMIs for barcode ranking is the best approach for most single-cell RNA-seq datasets. It generally performs well when the ambient RNA contamination is low to medium. For more contaminated datasets, it may be more robust to use the number of genes for barcode ranking. Generally, we suggest to test using the number of UMIs for barcode ranking first and visually inspect the plot If the threshold is not satisfactory, you may get better result using the number of genes for barcode ranking.
#'
#' \strong{Setting alpha}\cr
#' Breakpoints can only be found within the two vertical blacklines on the diagnostic plot. If the desired breakpoint is located outside of these lines, it is necessary to decrease alpha. If the desired breakpoint is far within the region, it is possible to increase alpha to improve convergence. If the selected alpha is too low, the function will automatically increment alpha untill it reaches alpha.max.
#' 
#' \strong{Setting boot}\cr
#' The function uses the root-mean-squared error to select the best segmentation model. The RMSE decreases with more breakpoints, therefore to choose a simple model that approximates the best model, the selected
#' @return A list or a data frame object that contains ranked barcodes and if chosen a threshold.
#' @export
#' @import segmented
#' @import zoo
#' @import Matrix

rank_barcodes = function(counts, type = "UMI", psi.min = 2, psi.max = 5, alpha = 0.001, alpha.max = 0.05, boot = 10, factor = 1.5, threshold = TRUE) {
  ## evaluate arguments
  # count matrix
  if(missing(counts)) {
    stop('No count matrix was provided', call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix","DelayedMatrix"))) { stop('Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix, DelayedMatrix', call. = FALSE) }
  }
  
  # type argument
  if(!any(type %in% c("Genes","gene","genes", "UMI", "umi","UMIS","umis","UMIs"))) stop('Incorrect input. Did you choose UMI or Genes?', call. = FALSE)
  
  # psi.min argument
  if(class(psi.min) != "numeric" | psi.min <= 0 | psi.min > psi.max) stop('psi.min needs to be a numeric greater than 0 and less than or equal to psi.max ', call. = FALSE)
  
  # psi.max argument
  if(class(psi.max) != "numeric" | psi.max <= 0 | psi.max < psi.min) stop('psi.max needs to be a numeric greater than 0 and greater than or equal to psi.min ', call. = FALSE)
  
  # threshold argument
  if(class(threshold) != "logical") stop('threshold needs to be a boolean (TRUE or FALSE)', call. = FALSE)
  
  ## convert the counts into dgCMatrix if its class() is not dgCMatrix
  if(class(counts) != "dgCMatrix") { counts = as(counts, "dgCMatrix") }
  
  ## get the feature type (allowing for spelling variants)
  feature_type <- "UMI"
  if (type %in% c("Genes","gene","genes")) { feature_type <- "Genes" }
  
  ## get barcode ranks using the appropriate type of feature
  switch(feature_type,
         UMI = {bcranks <- data.frame(counts = Matrix::colSums(counts))
         rownames(bcranks) <- colnames(counts)
         },
         Genes = {bcranks <- data.frame(counts = Matrix::colSums(counts > 0) )
         rownames(bcranks) <- colnames(counts)
         }
  )
  
  ## rank, filter and sort
  bcranks$rank <- rank(-bcranks$counts)
  bcranks <- bcranks[ bcranks$counts > 0,]
  bcranks <- bcranks[ order(-bcranks$counts, -bcranks$rank),]
  
  ## get unique counts, and log transform
  unique.counts <- bcranks[duplicated(bcranks$counts)==F,]
  unique.counts$counts <- log(unique.counts$counts)
  unique.counts$rank <- log(unique.counts$rank)
  
  ## Smooth the data using a moving average
  n <- ceiling(2*(nrow(unique.counts)^(1/3)))
  y <- zoo::rollmean(unique.counts$counts, k = n, align = "center")
  x <- zoo::rollmean(unique.counts$rank, k = n, align = "center")
  
  ## breakpoint analysis
  rmse <- as.data.frame(matrix(ncol=3, nrow = length(psi.min:psi.max)))
  models <- list()
  counter <- 1
  for (psi in psi.min:psi.max) {
    curr.alpha <- alpha
    model <- lm(y ~ x)
    out <- tryCatch(suppressWarnings(segmented::segmented(model, psi = seq(quantile(x, prob = curr.alpha),quantile(x, prob = (1-curr.alpha)), length.out = psi), control = segmented::seg.control(alpha = (curr.alpha-(curr.alpha/1000)), n.boot = boot))), error = function(e) e)
    if (class(out)[1] == "segmented") {
      rmse[counter,1] <- psi
      rmse[counter,2] <- sqrt(mean(out$residuals^2))
      rmse[counter,3] <- counter
      models[[counter]] <- out
      counter <- counter + 1
    } else if (any(grepl("psi values too close", out[[1]]))) {
      stop = 0
      while (stop == 0) {
        curr.alpha <- curr.alpha + alpha
        if (curr.alpha > alpha.max) { stop = 1 }
        out <- tryCatch(suppressWarnings(segmented::segmented(model, psi = seq(quantile(x, prob = curr.alpha),quantile(x, prob = (1-curr.alpha)), length.out = psi), control = segmented::seg.control(alpha = (curr.alpha-(curr.alpha/1000)), n.boot = boot))), error = function(e) e)
        if (class(out)[1] == "segmented") {
          rmse[counter,1] <- psi
          rmse[counter,2] <- sqrt(mean(out$residuals^2))
          rmse[counter,3] <- counter
          models[[counter]] <- out
          counter <- counter + 1
          stop = 1
        }
      }
    }
  }
  
  ## select the best model (within a factor of the smallest RMSE)
  rmse <- rmse[!is.na(rmse[,1]),]
  out <- models[[min(rmse[ rmse[,2] <= min(rmse[,2])*factor,3])]]
  
  ## select lower threshold
  slope <- segmented::slope(out)$x[,1]
  angles <- c()
  for (iter in 1:(length(slope) - 1)) { angles <- c(angles, atan((slope[iter]-slope[iter + 1]) / (1 + (slope[iter] * slope[iter + 1])))*(180/pi)) }
  best_bpt <- which.min(angles[ -1 ])+1
  lower_rank <- unique.counts[ which.min(abs(unique.counts$rank - out$psi[best_bpt,2])),2]
  lower <- unique.counts[ which.min(abs(unique.counts$rank  - out$psi[best_bpt,2])),1]
  
  ## finalize results depending on arguments
  if(threshold){
    output = list(ranks = bcranks, lower.threshold = exp(lower))
  } else {
    output = bcranks
  }
  
  ## return results
  return(output)
}

#' Identify Empty Droplets in Single-Cell RNA-seq Data
#'
#' @param mtx A sparse matrix format with features on the rows and barcodes on the columns (feature x barcodes)
#' @param fdr A numeric value specifying the false discovery rate (FDR) threshold for identifying empty droplets. Default is 0.001.
#'
#' @returns A list containing:
#'   - `rank`: Barcode ranking information.
#'   - `lower_threshold`: The lower threshold for barcode ranking.
#'   - `not_empty`: Vector of barcodes identified as non-empty droplets.
#' 
#' @examples identify_empty_droplets(mtx, fdr = 0.001)
empty_droplets <- function(mtx, fdr = 0.001, manual_thres = NULL, ...){
  # Get barcode ranks and threshold
  rank<- rank_barcodes(mtx, ...)
  
  if(!is.null(manual_thres)){
    # Find empty droplets using manual threshold
    empty <- DropletUtils::emptyDrops(mtx, 
                                      lower = manual_thres, 
                                      BPPARAM = BiocParallel::MulticoreParam())
    
    mtx_filter <-
      mtx[, BiocGenerics::colnames(mtx) %in%
            BiocGenerics::rownames(empty)[base::which(empty$FDR <= fdr)]]
    
    output <- base::list(rank = rank$ranks,
                         lower_threshold = data.frame(automatic = rank$lower.threshold,
                                                      manual = manual_thres),
                         not_empty = data.frame(not_empty = BiocGenerics::colnames(mtx_filter)))
  } else {
    # Find empty droplets using automatic threshold
    empty <- DropletUtils::emptyDrops(mtx, 
                                      lower = rank$lower.threshold, 
                                      BPPARAM = BiocParallel::MulticoreParam())
    mtx_filter <-
      mtx[, BiocGenerics::colnames(mtx) %in%
            BiocGenerics::rownames(empty)[base::which(empty$FDR <= fdr)]]
    
    output <- base::list(rank = rank$ranks,
                         lower_threshold = data.frame(automatic = rank$lower.threshold),
                         not_empty = data.frame(not_empty = BiocGenerics::colnames(mtx_filter)))
  }
  
  return(output)
}

identify_empty_droplets <- function(study_metadata, study_name, defaults = list(psi.max = 5, manual_thres = NULL), donor_settings = NULL, save_path, ...) {
  ## Identify empty drops
  results <-  study_metadata |> 
    ### Load mtx files to use for empty droplet detection
    dplyr::mutate(selected_mtx = purrr::pmap(list(cell_nuclei, name, sample), function(cell_nuclei, name, sample) {
      if (cell_nuclei == "cell") {
        read_mtx(path = paste0("/work/scRNAseq/", name, "/Preprocessed/", sample, "/Solo.out"), feature = "Gene") |> purrr::pluck(sample)
      } else if (cell_nuclei == "nuclei") {
        read_mtx(path = paste0("/work/scRNAseq/", name, "/Preprocessed/", sample, "/Solo.out"), feature = "GeneFull_Ex50pAS") |> purrr::pluck(sample)
      } else {
        stop("read of mtx failed")
      }
    }), 
      selected_mtx = purrr::map2(selected_mtx, ic_id, ~{colnames_mtx <- BiocGenerics::colnames(.x)
      BiocGenerics::colnames(.x) <- base::paste0(.y, "_", colnames_mtx)
      return(.x)
      }),
      ### Extract donor-specific settings 
      donor_params = purrr::map(ic_id, ~ {
        if (!is.null(donor_settings_droplet) && .x %in% base::names(donor_settings_droplet)) {
          # Merge donor-specific settings with defaults
          utils::modifyList(defaults,donor_settings_droplet[[.x]])
        } else {
          # Use defaults if no donor-specific settings exist
          defaults
        }
      }),
      ### Identify empty droplets
      barcodes = purrr::map2(selected_mtx, donor_params, function(mtx, params) {
        # Verify the input data type
        if (!inherits(mtx, "dgCMatrix")) {
          stop("Input is not a dgCMatrix")
        }
        # Call empty_droplets with verified inputs and explicit arguments
        empty_droplets(mtx, manual_thres = params$manual_thres, psi.max = params$psi.max)
      })) |> 
    dplyr::select(-selected_mtx)
  
  ## Save ranks
  ranks <- results |> dplyr::pull(barcodes) |> 
    purrr::map(purrr::pluck("rank")) |>
    purrr::set_names(results$ic_id)
  
  ### Save ranks as csv 
  purrr::imap(ranks, function(x, idx) {tibble::rownames_to_column(x, "barcode") |>  
      vroom::vroom_write(base::paste0(save_path, 
                                      study_name, 
                                      "_",
                                      idx,
                                      "_barcode_ranks.csv"),
                         delim = ",", col_names = TRUE)})
  
  ### Save threshold as csv
  threshold <- results |> dplyr::pull(barcodes) |> 
    purrr::map(purrr::pluck("lower_threshold")) |>
    purrr::set_names(results$ic_id)
  
  ### Save thresholds as csv 
  purrr::imap(threshold, function(x, idx) {x |>  
      vroom::vroom_write(base::paste0(save_path, 
                                      study_name, 
                                      "_",
                                      idx,
                                      "_lower_threshold.csv"),
                         delim = ",", col_names = TRUE)})
  
  ### Save non empty droplets as csv
  not_empty <- results |> dplyr::pull(barcodes) |> 
    purrr::map(purrr::pluck("not_empty")) |>
    purrr::set_names(results$ic_id)
  
  ### Save thresholds as csv 
  purrr::imap(not_empty, function(x, idx) {x |>  
      vroom::vroom_write(base::paste0(save_path, 
                                      study_name, 
                                      "_",
                                      idx,
                                      "_not_empty.csv"),
                         delim = ",", col_names = TRUE)})
  return(NULL)
}


#' Compute Quality Metrics for Single-Cell RNA-seq Data
#'
#' @param mtx A sparse matrix format with features on the rows and barcodes on the columns (feature x barcodes)
#' @param mitogenes A vector of mitochondrial gene names
#' @param ribogenes A vector of ribosomal gene names
#' @param pcgenes A vector of protein-coding gene names
#' @param mtx_gene A sparse matrix of exon counts (feature x barcodes)
#' @param mtx_genefull A sparse matrix of exon+intron counts (feature x barcodes)
#'
#' @returns A dataframe with quality control metrics
#' @examples quality_metrics(mtx, mitogenes, ribogenes, pcgenes, mtx_gene, mtx_genefull)
#' 
quality_metrics <- function(mtx, mitogenes, ribogenes, pcgenes, mtx_gene, mtx_genefull){
  
  ## Contrast 
  # Match barcodes between the input counts and the mtx_gene matrix
  mtx_gene_sub <- mtx_gene[, BiocGenerics::colnames(mtx_gene) %in% BiocGenerics::colnames(mtx)]
  mtx_gene_sub <- mtx_gene_sub[, base::match(BiocGenerics::colnames(mtx), BiocGenerics::colnames(mtx_gene_sub))]
  # Match barcodes between the input counts and the mtx_gene matrix
  mtx_genefull_sub <- mtx_genefull[, BiocGenerics::colnames(mtx_genefull) %in% BiocGenerics::colnames(mtx)]
  mtx_genefull_sub <- mtx_genefull_sub[, base::match(BiocGenerics::colnames(mtx), BiocGenerics::colnames(mtx_genefull_sub))]
  
  # Internalize matrix
  metrics <- BiocGenerics::as.data.frame(base::matrix(ncol=8, 
                                                      nrow=ncol(mtx)))
  # Set colnames
  BiocGenerics::colnames(metrics) <- c("barcode",
                                       "logUMIs",
                                       "logFeatures",
                                       "mitochondrial_fraction",
                                       "ribosomal_fraction",
                                       "coding_fraction", 
                                       "contrast_fraction", 
                                       "complexity")
  
  metrics[,1] <- BiocGenerics::colnames(mtx)
  metrics[,2] <- Matrix::colSums(mtx)
  metrics[,3] <- Matrix::colSums(mtx > 0)
  metrics[,4] <- Matrix::colSums(mtx[ BiocGenerics::rownames(mtx) %in% mitogenes,]) / metrics[,2]
  metrics[,5] <- Matrix::colSums(mtx[ BiocGenerics::rownames(mtx) %in% ribogenes,]) / metrics[,2]
  metrics[,6] <- Matrix::colSums(mtx[ BiocGenerics::rownames(mtx) %in% pcgenes,]) / metrics[,2]
  metrics[,7] <- Matrix::colSums(mtx_gene_sub) / Matrix::colSums(mtx_genefull_sub)
  metrics[,8] <- metrics[,3] / metrics[,2]
  metrics[,2] <- base::log(metrics[,2])
  metrics[,3] <- base::log(metrics[,3])
  
  return(metrics)
}

quality_metrics_per_sample <- function(path, study_metadata, genes, study_name, not_empty) {
  
  # Extract unique library prep
  library_prep <- base::unique(study_metadata$library_prep)
  base::message("library_prep: '", library_prep, "'")
  
  platform <- base::unique(study_metadata$platform)
  
  if (length(platform) != 1) stop("Multiple library preps found in study. Please split them first.")
  
  # Define expected order of samples
  expected_ids <- study_metadata$ic_id
  
  # Validate input length consistency
  if (base::length(study_metadata$cell_nuclei) != base::length(study_metadata$ic_id)) {
    stop("Mismatch in lengths of `cell_nuclei` and `ic_id`.")
  }
  
  ## Load all mtx files 
  mtx_gene_list <- read_mtx(path, "Gene")
  mtx_gene_list <- sample_to_ic_id(mtx_gene_list, study_metadata)
  mtx_gene_list <- prefix_colnames(mtx_gene_list)
  
  mtx_genefull_list <- read_mtx(path, "GeneFull")
  mtx_genefull_list <- sample_to_ic_id(mtx_genefull_list, study_metadata)
  mtx_genefull_list <- prefix_colnames(mtx_genefull_list)
  
  # Make sure files are in the expected order
  mtx_gene_list <- mtx_gene_list[base::match(expected_ids, base::names(mtx_gene_list))]
  mtx_genefull_list <- mtx_genefull_list[base::match(expected_ids, base::names(mtx_genefull_list))]
  
  ### Define count methods which will be used 
  rna_count <- "exon"
  
  # Load GeneFull_Ex50pAS only if any sample is nuclei
  mtx_genefull_50_list <- NULL
  if ("nuclei" %in% study_metadata$cell_nuclei) {
    base::message("Nuclei samples, loading GeneFull_Ex50pAS")
    
    mtx_genefull_50_list <- read_mtx(path, "GeneFull_Ex50pAS")
    mtx_genefull_50_list <- sample_to_ic_id(mtx_genefull_50_list, study_metadata)
    mtx_genefull_50_list <- prefix_colnames(mtx_genefull_50_list)
    
    # Make sure files are in the expected order
    mtx_genefull_50_list <- mtx_genefull_50_list[base::match(expected_ids, base::names(mtx_genefull_50_list))]
    
    
    ### Define count methods which will be used 
    rna_count <- "exon_intron"
  }
  
  # Check that all required ic_ids are present in the loaded lists
  missing_gene_ids <- BiocGenerics::setdiff(expected_ids, names(mtx_gene_list))
  if (base::length(missing_gene_ids) > 0) {
    stop(base::paste("The following ic_ids are missing in `mtx_gene_list`:", base::paste(missing_gene_ids, collapse = ", ")))
  }
  if (!is.null(mtx_genefull_50_list)) {
    missing_50_ids <- BiocGenerics::setdiff(expected_ids, base::names(mtx_genefull_50_list))
    if (base::length(missing_50_ids) > 0) {
      stop(base::paste("The following ic_ids are missing in `mtx_genefull_50_list`:", base::paste(missing_50_ids, collapse = ", ")))
    }
  }
  
  # Check if it's droplet, plate-based or plate-based with barcode
  if (base::as.character(platform) == "droplet") {
    base::message("Method is droplet based")
    
    # Get not empty barcodes
    not_empty_barcodes <- not_empty[grep(study_name, names(not_empty))] |> 
      dplyr::bind_rows() |> dplyr::pull("not_empty")
    
    ## Droplet quality control ----
    results <- study_metadata %>%
      dplyr::mutate(
        ### Select main mtx file to use 
        selected_mtx = purrr::map2(cell_nuclei, ic_id, ~ {
          if (.x == "cell" && .y %in% names(mtx_gene_list)) {
            base::message(paste("Selected mtx_gene for ic_id:", .y))
            mtx_gene_list[[.y]]
          } else if (.y %in% base::names(mtx_genefull_50_list)) {
            base::message(base::paste("Selected mtx_genefull_50 for ic_id:", .y))
            mtx_genefull_50_list[[.y]]
          } else {
            stop(base::message(base::paste("No valid matrix found for ic_id:", .y)))
          }
        }),
        rna_count = rna_count,
        ### Extract non empty droplets
        selected_mtx = purrr::map(selected_mtx, ~.x[, BiocGenerics::colnames(.x) %in% not_empty_barcodes]),
        ### Quality control metrics 
        quality_met = purrr::pmap(list(selected_mtx, mtx_gene_list, mtx_genefull_list),
                                  ~ quality_metrics(
                                    mtx = ..1,
                                    mitogenes = genes[["mito_genes"]],
                                    ribogenes = genes[["ribo_genes"]],
                                    pcgenes = genes[["protein_genes"]],
                                    mtx_gene = ..2,
                                    mtx_genefull = ..3
                                  ))) |> 
      dplyr::select(-selected_mtx) |> 
      dplyr::relocate(rna_count, quality_met, .after = library_prep)
    
    ### Save quality control metrics 
    results <- results |> 
      tidyr::unnest(quality_met)
    
    
    vroom::vroom_write(results, base::paste0(here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics/"),
                                             study_name,
                                             "_quality_metrics.csv"),
                       delim = ",", col_names = TRUE)
    
    
  } else if (base::as.character(platform) == "plate_barcode") {
    base::message("Method is plate based with barcode")
    
    ## Plate-based barcode ----
    results <- study_metadata %>%
      dplyr::mutate(
        ### Select main mtx file to use 
        selected_mtx = purrr::map2(cell_nuclei, ic_id, ~ {
          if (.x == "cell" && .y %in% names(mtx_gene_list)) {
            base::message(paste("Selected mtx_gene for ic_id:", .y))
            mtx_gene_list[[.y]]
          } else if (.y %in% base::names(mtx_genefull_50_list)) {
            base::message(base::paste("Selected mtx_genefull_50 for ic_id:", .y))
            mtx_genefull_50_list[[.y]]
          } else {
            stop(base::message(base::paste("No valid matrix found for ic_id:", .y)))
          }
        }),
        rna_count = rna_count,
        ### Quality control metrics 
        quality_met = purrr::pmap(list(selected_mtx, mtx_gene_list, mtx_genefull_list),
                                  ~ quality_metrics(
                                    mtx = ..1,
                                    mitogenes = genes[["mito_genes"]],
                                    ribogenes = genes[["ribo_genes"]],
                                    pcgenes = genes[["protein_genes"]],
                                    mtx_gene = ..2,
                                    mtx_genefull = ..3
                                  ))) |> 
      dplyr::select(-selected_mtx) |> 
      dplyr::relocate(rna_count, quality_met, .after = library_prep)
    
    ### Save quality control metrics 
    results <- results |> 
      tidyr::unnest(quality_met)
    
    vroom::vroom_write(results, base::paste0(here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics/"),
                                             study_name,
                                             "_quality_metrics.csv"),
                       delim = ",", col_names = TRUE)
    
  } else if (base::as.character(platform) == "plate") {
    base::message("Method is plate based")
    
    ## Plate-based ----
    
    ## Combine the matrices column-wise into a single large sparse matrix 
    mtx_gene <- base::do.call(base::cbind, mtx_gene_list)
    mtx_genefull <- base::do.call(base::cbind, mtx_genefull_list)
    
    if (!is.null(mtx_genefull_50_list)) {
      mtx_genefull_50 <- base::do.call(base::cbind, mtx_genefull_50_list) 
    }
    
    study_metadata <- study_metadata |> 
      dplyr::mutate(barcode = ic_id,
                    rna_count = rna_count)
    
    # Initialize empty lists for quality metrics
    quality_met_nuclei <- NULL
    quality_met_cells <- NULL
    
    # Separate cells and nuclei
    nuclei_metadata <- study_metadata |> dplyr::filter(cell_nuclei == "nuclei")
    cell_metadata <- study_metadata |> dplyr::filter(cell_nuclei == "cell")
    
    # Calculate QC metrics separately
    if (nrow(nuclei_metadata) > 0) {
      base::message("Processing QC for nuclei...")
      quality_met_nuclei <- quality_metrics(
        mtx = mtx_genefull_50,
        mitogenes = genes[["mito_genes"]],
        ribogenes = genes[["ribo_genes"]],
        pcgenes = genes[["protein_genes"]],
        mtx_gene = mtx_gene,
        mtx_genefull = mtx_genefull
      )
    }
    
    if (nrow(cell_metadata) > 0) {
      base::message("Processing QC for cells...")
      quality_met_cells <- quality_metrics(
        mtx = mtx_gene,
        mitogenes = genes[["mito_genes"]],
        ribogenes = genes[["ribo_genes"]],
        pcgenes = genes[["protein_genes"]],
        mtx_gene = mtx_gene,
        mtx_genefull = mtx_genefull
      )
    }
    
    # Combine results
    # Check if both dataframes exist, and merge them only if both are available
    if (!is.null(quality_met_nuclei) & !is.null(quality_met_cells)) {
      quality_met <- dplyr::bind_rows(quality_met_nuclei, quality_met_cells)
    } else if (!is.null(quality_met_nuclei)) {
      quality_met <- quality_met_nuclei
    } else if (!is.null(quality_met_cells)) {
      quality_met <- quality_met_cells
    } else {
      stop("No quality metrics calculated for either nuclei or cells.")
    }
    
    ## Merge results
    results <- dplyr::full_join(x = study_metadata, y = quality_met, by = "barcode") |> 
      dplyr::relocate(rna_count, barcode,
                      logUMIs,
                      logFeatures,
                      mitochondrial_fraction,
                      ribosomal_fraction,
                      coding_fraction, 
                      contrast_fraction, 
                      complexity, .after = library_prep)
    
    ## Save CSV  
    vroom::vroom_write(
      results,
      base::paste0(here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics/"),
                   study_name,
                   "_quality_metrics.csv"),
      delim = ",",
      col_names = TRUE)
    
  } else {
    stop("Unknown library preparation type")
  }
  
  return(NULL)
}


read_gtf_parallel <- function(file, max_cores = NULL, chunk_size = 10, max_size = 2 * 1024^3) {
  # Load necessary libraries
  
  # Set up the parallel backend
  if (is.null(max_cores)) {
    plan(multisession) # Use all available cores by default
  } else {
    plan(multisession, workers = max_cores) # Limit cores if specified
  }
  
  # Set the maximum size for globals in parallel processing
  options(future.globals.maxSize = max_size)
  
  # Column names for GTF
  cnames <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  
  # Read the GTF file
  messy <- read_tsv(file, col_names = cnames, comment = "#")
  
  # Extract unique attribute names
  att_names <- str_split(messy$attribute[1:100], '"; ') |> 
    unlist() |> 
    trimws() |> 
    trimws(whitespace = ";") |> 
    sub(pattern = " .*$", replacement = "") |> 
    unique()
  
  att_names <- att_names[att_names != ""]
  
  # Split attribute names into chunks
  chunks <- split(att_names, ceiling(seq_along(att_names) / chunk_size))
  
  # Process each chunk in parallel
  results <- map_dfc(chunks, function(chunk) {
    future_map_dfc(chunk, function(att) {
      extracted <- str_extract(messy$attribute, sprintf('";\\s+%1$s[^;]+|^%1$s[^;]+;[^"]+"', att)) |>
        trimws(whitespace = '["; ]+', which = 'left') |>
        str_extract('(?<=")[^"]+(?=")')
      tibble("{att}" := extracted)
    })
  })
  
  # Combine results with the original data (removing the attribute column)
  messy <- bind_cols(messy, results) |> select(-attribute)
  
  return(messy)
}



# plotting functions ------------------------------------------------------
my_theme <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(strip.background = element_blank(),
                   plot.title = element_text(size = 8, hjust = 0.5),
                   axis.title = element_text(size = 7),
                   axis.text = element_text(size = 6, color = "black"),
                   axis.ticks = element_line(color = "black"),
                   legend.text = element_text(size = 4), 
                   legend.title = element_text(size = 4),
                   panel.background = element_rect(fill='transparent'),
                   plot.background = element_rect(fill='transparent', color=NA), 
                   legend.background = element_rect(fill='transparent', color = NA), 
                   legend.box.background = element_rect(fill='transparent', color = NA))
}

plot_rank <- function(rank, threshold, title, non_empty){
  p <- rank |> 
    dplyr::mutate(not_empty = dplyr::case_when(barcode %in% non_empty$not_empty ~ "Not empty",
                                               .default = "Empty"),
                  not_empty = factor(not_empty, levels = c("Not empty", "Empty"))) |> 
    ggplot2::ggplot(aes(x = rank, y = counts, color = not_empty)) +
    ggrastr::geom_point_rast(size = 0.1, shape = 20, raster.dpi = 500) +
    ggplot2::scale_color_manual(values = c("Not empty" = "black", "Empty" = "#D3D3D3")) +
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggplot2::labs(x = "Log Rank",
                  y = "Log UMIs", 
                  color = "Droplet",
                  linetype = "Lower threshold",
                  title = title) +
    ggplot2::geom_hline(aes(yintercept = threshold$automatic, linetype = "Automatic"), color = "red") +
    ggplot2::scale_linetype_manual(values = c("Automatic" = 3, "Manual" = 5)) +
    my_theme() +
    ggplot2::theme(aspect.ratio=1,
                   legend.position = "inside",
                   legend.position.inside =  c(0.2, 0.2),
                   legend.key.size = unit(0.1, "inch"),
                   legend.spacing.y = unit(0, "inch"),
                   legend.margin = margin(0, 0, 0, 0),
                   legend.text = element_text(margin = margin(l = 0)),
                   legend.title = element_blank())
  
  if (!is.null(threshold$manual)) {
   p <- p +  ggplot2::geom_hline(aes(yintercept = threshold$manual, linetype = "Manual"), color = "blue")
  }
  
  return(p)
}

plot_rank_line <- function(rank, title, failed){
  p <- rank |> 
    dplyr::mutate(Passed = dplyr::case_when(name %in% failed_samples ~ "Failed",
                                            .default = "Passed")) |> 
    ggplot2::ggplot(aes(x = rank, y = counts, group = name, color = Passed)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values= c("Passed" = "black", "Failed" = "red")) +
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggplot2::labs(x = "Log Rank",
                  y = "Log UMIs", 
                  title = title) +
    my_theme() +
    ggplot2::theme(aspect.ratio=1,
                   legend.position = "inside",
                   legend.position.inside =  c(0.2, 0.2),
                   legend.key.size = unit(0.1, "inch"),
                   legend.spacing.y = unit(0, "inch"),
                   legend.margin = margin(0, 0, 0, 0),
                   legend.text = element_text(margin = margin(l = 0)),
                   legend.title = element_blank())
  
  return(p)
}

plot_hist_qc <- function(.x, .y){
  ggplot2::ggplot(data = tibble::tibble(value = .x), ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 100, fill = "black") +
    ggplot2::labs(title = dplyr::case_when(
      .y == "logUMIs" ~ "Number of unique transcripts",
      .y == "logFeatures" ~ "Number of detected features",
      .y == "mitochondrial_fraction" ~ "Mitochondrial fraction",
      .y == "ribosomal_fraction" ~ "Ribosomal fraction",
      .y == "coding_fraction" ~ "Protein coding fraction",
      .y == "contrast_fraction" ~ "Contrast fraction",
      .y == "complexity" ~ "Library complexity",
      TRUE ~ "value"), 
      x = dplyr::case_when(
        .y == "logUMIs" ~ "LogUMIs",
        .y == "logFeatures" ~ "LogFeatures",
        .y == "mitochondrial_fraction" ~ "% of counts in mitochondrial genes",
        .y == "ribosomal_fraction" ~ "% of counts in ribosomal genes",
        .y == "coding_fraction" ~ "% of counts in protein coding genes",
        .y == "contrast_fraction" ~ "exon / exon+intron counts",
        .y == "complexity" ~ "LogFeatures/LogUMI",
        TRUE ~ "value"), 
      y = "Frequency") +
    my_theme() +
    ggplot2::theme(aspect.ratio=1)
}

plot_hist_star <- function(.x, .y){
  ggplot2::ggplot(data = tibble::tibble(value = .x), ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 100, fill = "black") +
    ggplot2::labs(title = dplyr::case_when(
      .y == "logUMIs" ~ "Number of unique transcripts",
      .y == "logFeatures" ~ "Number of detected features",
      .y == "mitochondrial_fraction" ~ "Mitochondrial fraction",
      .y == "ribosomal_fraction" ~ "Ribosomal fraction",
      .y == "coding_fraction" ~ "Protein coding fraction",
      .y == "contrast_fraction" ~ "Contrast fraction",
      .y == "complexity" ~ "Library complexity",
      TRUE ~ "value"), 
      x = dplyr::case_when(
        .y == "logUMIs" ~ "LogUMIs",
        .y == "logFeatures" ~ "LogFeatures",
        .y == "mitochondrial_fraction" ~ "% of counts in mitochondrial genes",
        .y == "ribosomal_fraction" ~ "% of counts in ribosomal genes",
        .y == "coding_fraction" ~ "% of counts in protein coding genes",
        .y == "contrast_fraction" ~ "exon / exon+intron counts",
        .y == "complexity" ~ "LogFeatures/LogUMI",
        TRUE ~ "value"), 
      y = "Frequency") +
    my_theme() +
    ggplot2::theme(aspect.ratio=1)
}
