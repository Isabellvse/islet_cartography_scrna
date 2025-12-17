#' Create directories if they do not exist
#'
#' This function checks if each directory in the provided list exists. If any directory
#' does not exist, it creates the directory (along with any necessary parent directories).
#' It prints a message for each directory indicating whether it was created or already exists.
#'
#' @param output_dir A character vector or list of directory paths to be checked and created.
#'
#' @return NULL This function performs a side-effect (creating directories) and does not return anything.
create_directories <- function(output_dir) {
  purrr::walk(output_dir, ~ {
    if (!dir.exists(.x)) {
      dir.create(.x, recursive = TRUE) # create the directory if it doesn't exist
      print(paste0(.x, " has been created!"))
    } else {
      print(paste0(.x, " already exists!"))
    }
  })
}


# Loading  ----------------------------------------------------------------
load_data_with_classes <- function(data_file = "data.csv", class_file = "column_classes.csv") {
  # Read the column classes
  df_classes <- vroom::vroom(class_file, delim = ",", col_types = cols(
    column = vroom::col_character(),
    class = vroom::col_character()
  ))

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
        vroom::col_guess() # Default to guessing for unknown types
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
load_star_quality <- function(path, i) {
  df <- path |>
    vroom::vroom() |>
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::contains("%"),
        .fns = ~ stringr::str_remove(.x, pattern = "%") |>
          base::as.numeric()
      ),
      dplyr::across(
        .cols = c("Deletion_rate_per_base", "Insertion_rate_per_base"),
        .fns = ~ stringr::str_remove(.x, pattern = "%") |>
          base::as.numeric()
      ),
      name = i
    ) |>
    dplyr::rename(
      "Deletion_rate_per_base_%" = "Deletion_rate_per_base",
      "Insertion_rate_per_base_%" = "Insertion_rate_per_base"
    ) |>
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
read_mtx <- function(path, feature) {
  # Get list of all relevant files
  base_path <- base::paste0(path, "/", feature, "/", "raw")

  base::message(base::paste0("Reading mtx files from: ", base_path))
  # Get file paths
  files <- base::list.files(
    path = base::Sys.glob(base_path),
    pattern = "matrix.mtx|barcodes.tsv|features.tsv",
    recursive = TRUE,
    full.names = TRUE
  )

  # Extract sample names from file paths
  file_info <- base::data.frame(
    path = files,
    sample = stringr::str_extract(files, "(?<=Preprocessed/)[^/]+"),
    matrix = stringr::str_extract(files, "barcodes.tsv|features.tsv|matrix.mtx")
  )

  # Reorder dataframe to ensure correct matching of paths
  file_group <- file_info |>
    tidyr::pivot_wider(names_from = matrix, values_from = path) |>
    dplyr::rename(
      mtx = `matrix.mtx`,
      barcodes = `barcodes.tsv`,
      features = `features.tsv`
    ) |>
    tibble::column_to_rownames("sample")

  mtx_list <- purrr::pmap(file_group, function(mtx, barcodes, features) {
    Seurat::ReadMtx(mtx = mtx, cells = barcodes, features = features, feature.column = 1)
  }) |>
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

rank_barcodes <- function(counts, type = "UMI", psi.min = 2, psi.max = 5, alpha = 0.001, alpha.max = 0.05, boot = 10, factor = 1.5, threshold = TRUE) {
  ## evaluate arguments
  # count matrix
  if (missing(counts)) {
    stop("No count matrix was provided", call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix", "matrix", "dgCMatrix", "DelayedMatrix"))) {
      stop("Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix, DelayedMatrix", call. = FALSE)
    }
  }

  # type argument
  if (!any(type %in% c("Genes", "gene", "genes", "UMI", "umi", "UMIS", "umis", "UMIs"))) stop("Incorrect input. Did you choose UMI or Genes?", call. = FALSE)

  # psi.min argument
  if (class(psi.min) != "numeric" | psi.min <= 0 | psi.min > psi.max) stop("psi.min needs to be a numeric greater than 0 and less than or equal to psi.max ", call. = FALSE)

  # psi.max argument
  if (class(psi.max) != "numeric" | psi.max <= 0 | psi.max < psi.min) stop("psi.max needs to be a numeric greater than 0 and greater than or equal to psi.min ", call. = FALSE)

  # threshold argument
  if (class(threshold) != "logical") stop("threshold needs to be a boolean (TRUE or FALSE)", call. = FALSE)

  ## convert the counts into dgCMatrix if its class() is not dgCMatrix
  if (class(counts) != "dgCMatrix") {
    counts <- as(counts, "dgCMatrix")
  }

  ## get the feature type (allowing for spelling variants)
  feature_type <- "UMI"
  if (type %in% c("Genes", "gene", "genes")) {
    feature_type <- "Genes"
  }

  ## get barcode ranks using the appropriate type of feature
  switch(feature_type,
    UMI = {
      bcranks <- data.frame(counts = Matrix::colSums(counts))
      rownames(bcranks) <- colnames(counts)
    },
    Genes = {
      bcranks <- data.frame(counts = Matrix::colSums(counts > 0))
      rownames(bcranks) <- colnames(counts)
    }
  )

  ## rank, filter and sort
  bcranks$rank <- rank(-bcranks$counts)
  bcranks <- bcranks[bcranks$counts > 0, ]
  bcranks <- bcranks[order(-bcranks$counts, -bcranks$rank), ]

  ## get unique counts, and log transform
  unique.counts <- bcranks[duplicated(bcranks$counts) == F, ]
  unique.counts$counts <- log(unique.counts$counts)
  unique.counts$rank <- log(unique.counts$rank)

  ## Smooth the data using a moving average
  n <- ceiling(2 * (nrow(unique.counts)^(1 / 3)))
  y <- zoo::rollmean(unique.counts$counts, k = n, align = "center")
  x <- zoo::rollmean(unique.counts$rank, k = n, align = "center")

  ## breakpoint analysis
  rmse <- as.data.frame(matrix(ncol = 3, nrow = length(psi.min:psi.max)))
  models <- list()
  counter <- 1
  for (psi in psi.min:psi.max) {
    curr.alpha <- alpha
    model <- lm(y ~ x)
    out <- tryCatch(suppressWarnings(segmented::segmented(model, psi = seq(quantile(x, prob = curr.alpha), quantile(x, prob = (1 - curr.alpha)), length.out = psi), control = segmented::seg.control(alpha = (curr.alpha - (curr.alpha / 1000)), n.boot = boot))), error = function(e) e)
    if (class(out)[1] == "segmented") {
      rmse[counter, 1] <- psi
      rmse[counter, 2] <- sqrt(mean(out$residuals^2))
      rmse[counter, 3] <- counter
      models[[counter]] <- out
      counter <- counter + 1
    } else if (any(grepl("psi values too close", out[[1]]))) {
      stop <- 0
      while (stop == 0) {
        curr.alpha <- curr.alpha + alpha
        if (curr.alpha > alpha.max) {
          stop <- 1
        }
        out <- tryCatch(suppressWarnings(segmented::segmented(model, psi = seq(quantile(x, prob = curr.alpha), quantile(x, prob = (1 - curr.alpha)), length.out = psi), control = segmented::seg.control(alpha = (curr.alpha - (curr.alpha / 1000)), n.boot = boot))), error = function(e) e)
        if (class(out)[1] == "segmented") {
          rmse[counter, 1] <- psi
          rmse[counter, 2] <- sqrt(mean(out$residuals^2))
          rmse[counter, 3] <- counter
          models[[counter]] <- out
          counter <- counter + 1
          stop <- 1
        }
      }
    }
  }

  ## select the best model (within a factor of the smallest RMSE)
  rmse <- rmse[!is.na(rmse[, 1]), ]
  out <- models[[min(rmse[rmse[, 2] <= min(rmse[, 2]) * factor, 3])]]

  ## select lower threshold
  slope <- segmented::slope(out)$x[, 1]
  angles <- c()
  for (iter in 1:(length(slope) - 1)) {
    angles <- c(angles, atan((slope[iter] - slope[iter + 1]) / (1 + (slope[iter] * slope[iter + 1]))) * (180 / pi))
  }
  best_bpt <- which.min(angles[-1]) + 1
  lower_rank <- unique.counts[which.min(abs(unique.counts$rank - out$psi[best_bpt, 2])), 2]
  lower <- unique.counts[which.min(abs(unique.counts$rank - out$psi[best_bpt, 2])), 1]

  ## finalize results depending on arguments
  if (threshold) {
    output <- list(ranks = bcranks, lower.threshold = exp(lower))
  } else {
    output <- bcranks
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
empty_droplets <- function(mtx, fdr = 0.001, manual_thres = NULL, ...) {
  # Get barcode ranks and threshold
  rank <- rank_barcodes(mtx, ...)

  if (!is.null(manual_thres)) {
    # Find empty droplets using manual threshold
    empty <- DropletUtils::emptyDrops(mtx,
      lower = manual_thres,
      BPPARAM = BiocParallel::MulticoreParam()
    )

    mtx_filter <-
      mtx[, BiocGenerics::colnames(mtx) %in%
        BiocGenerics::rownames(empty)[base::which(empty$FDR <= fdr)]]

    output <- base::list(
      rank = rank$ranks,
      lower_threshold = data.frame(
        automatic = rank$lower.threshold,
        manual = manual_thres
      ),
      not_empty = data.frame(not_empty = BiocGenerics::colnames(mtx_filter))
    )
  } else {
    # Find empty droplets using automatic threshold
    empty <- DropletUtils::emptyDrops(mtx,
      lower = rank$lower.threshold,
      BPPARAM = BiocParallel::MulticoreParam()
    )
    mtx_filter <-
      mtx[, BiocGenerics::colnames(mtx) %in%
        BiocGenerics::rownames(empty)[base::which(empty$FDR <= fdr)]]

    output <- base::list(
      rank = rank$ranks,
      lower_threshold = data.frame(automatic = rank$lower.threshold),
      not_empty = data.frame(not_empty = BiocGenerics::colnames(mtx_filter))
    )
  }

  return(output)
}

identify_empty_droplets <- function(study_metadata, study_name, defaults = list(psi.max = 5, manual_thres = NULL), donor_settings = NULL, save_path, ...) {
  ## Identify empty drops
  results <- study_metadata |>
    ### Load mtx files to use for empty droplet detection
    dplyr::mutate(
      selected_mtx = purrr::pmap(list(cell_nuclei, name, sample), function(cell_nuclei, name, sample) {
        if (cell_nuclei == "cell") {
          read_mtx(path = paste0("/work/scRNAseq/", name, "/Preprocessed/", sample, "/Solo.out"), feature = "Gene") |> purrr::pluck(sample)
        } else if (cell_nuclei == "nuclei") {
          read_mtx(path = paste0("/work/scRNAseq/", name, "/Preprocessed/", sample, "/Solo.out"), feature = "GeneFull_Ex50pAS") |> purrr::pluck(sample)
        } else {
          stop("read of mtx failed")
        }
      }),
      selected_mtx = purrr::map2(selected_mtx, ic_id, ~ {
        colnames_mtx <- BiocGenerics::colnames(.x)
        BiocGenerics::colnames(.x) <- base::paste0(.y, "_", colnames_mtx)
        return(.x)
      }),
      ### Extract donor-specific settings
      donor_params = purrr::map(ic_id, ~ {
        if (!is.null(donor_settings_droplet) && .x %in% base::names(donor_settings_droplet)) {
          # Merge donor-specific settings with defaults
          utils::modifyList(defaults, donor_settings_droplet[[.x]])
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
      })
    ) |>
    dplyr::select(-selected_mtx)

  ## Save ranks
  ranks <- results |>
    dplyr::pull(barcodes) |>
    purrr::map(purrr::pluck("rank")) |>
    purrr::set_names(results$ic_id)

  ### Save ranks as csv
  purrr::imap(ranks, function(x, idx) {
    tibble::rownames_to_column(x, "barcode") |>
      vroom::vroom_write(
        base::paste0(
          save_path,
          study_name,
          "_",
          idx,
          "_barcode_ranks.csv"
        ),
        delim = ",", col_names = TRUE
      )
  })

  ### Save threshold as csv
  threshold <- results |>
    dplyr::pull(barcodes) |>
    purrr::map(purrr::pluck("lower_threshold")) |>
    purrr::set_names(results$ic_id)

  ### Save thresholds as csv
  purrr::imap(threshold, function(x, idx) {
    x |>
      vroom::vroom_write(
        base::paste0(
          save_path,
          study_name,
          "_",
          idx,
          "_lower_threshold.csv"
        ),
        delim = ",", col_names = TRUE
      )
  })

  ### Save non empty droplets as csv
  not_empty <- results |>
    dplyr::pull(barcodes) |>
    purrr::map(purrr::pluck("not_empty")) |>
    purrr::set_names(results$ic_id)

  ### Save thresholds as csv
  purrr::imap(not_empty, function(x, idx) {
    x |>
      vroom::vroom_write(
        base::paste0(
          save_path,
          study_name,
          "_",
          idx,
          "_not_empty.csv"
        ),
        delim = ",", col_names = TRUE
      )
  })
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
quality_metrics <- function(mtx, mitogenes, ribogenes, pcgenes, mtx_gene, mtx_genefull) {
  ## Contrast
  # Match barcodes between the input counts and the mtx_gene matrix
  mtx_gene_sub <- mtx_gene[, BiocGenerics::colnames(mtx_gene) %in% BiocGenerics::colnames(mtx)]
  mtx_gene_sub <- mtx_gene_sub[, base::match(BiocGenerics::colnames(mtx), BiocGenerics::colnames(mtx_gene_sub))]
  # Match barcodes between the input counts and the mtx_gene matrix
  mtx_genefull_sub <- mtx_genefull[, BiocGenerics::colnames(mtx_genefull) %in% BiocGenerics::colnames(mtx)]
  mtx_genefull_sub <- mtx_genefull_sub[, base::match(BiocGenerics::colnames(mtx), BiocGenerics::colnames(mtx_genefull_sub))]

  # Internalize matrix
  metrics <- BiocGenerics::as.data.frame(base::matrix(
    ncol = 10,
    nrow = ncol(mtx)
  ))
  # Set colnames
  BiocGenerics::colnames(metrics) <- c(
    "barcode",
    "nUMIs",
    "nFeatures",
    "lnUMIs",
    "lnFeatures",
    "mitochondrial_fraction",
    "ribosomal_fraction",
    "coding_fraction",
    "contrast_fraction",
    "complexity"
  )

  metrics[, 1] <- BiocGenerics::colnames(mtx)
  metrics[, 2] <- Matrix::colSums(mtx)
  metrics[, 3] <- Matrix::colSums(mtx > 0)
  metrics[, 4] <- base::log(metrics[, 2])
  metrics[, 5] <- base::log(metrics[, 3])
  metrics[, 6] <- Matrix::colSums(mtx[BiocGenerics::rownames(mtx) %in% mitogenes, ]) / metrics[, 2]
  metrics[, 7] <- Matrix::colSums(mtx[BiocGenerics::rownames(mtx) %in% ribogenes, ]) / metrics[, 2]
  metrics[, 8] <- Matrix::colSums(mtx[BiocGenerics::rownames(mtx) %in% pcgenes, ]) / metrics[, 2]
  metrics[, 9] <- Matrix::colSums(mtx_gene_sub) / Matrix::colSums(mtx_genefull_sub)
  metrics[, 10] <- metrics[, 3] / metrics[, 2]

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
      dplyr::bind_rows() |>
      dplyr::pull("not_empty")

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
        selected_mtx = purrr::map(selected_mtx, ~ .x[, BiocGenerics::colnames(.x) %in% not_empty_barcodes]),
        ### Quality control metrics
        quality_met = purrr::pmap(
          list(selected_mtx, mtx_gene_list, mtx_genefull_list),
          ~ quality_metrics(
            mtx = ..1,
            mitogenes = genes[["mito_genes"]],
            ribogenes = genes[["ribo_genes"]],
            pcgenes = genes[["protein_genes"]],
            mtx_gene = ..2,
            mtx_genefull = ..3
          )
        )
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate("Unmapped_reads_%" = sum(dplyr::c_across(tidyselect::starts_with("%_of_reads_unmapped_")), na.rm = TRUE)) %>%
      dplyr::ungroup() |>
      dplyr::select(-selected_mtx) |>
      dplyr::relocate(rna_count, quality_met, .after = library_prep)

    ### Save quality control metrics
    results <- results |>
      tidyr::unnest(quality_met)


    vroom::vroom_write(results, base::paste0(
      here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics/"),
      study_name,
      "_quality_metrics.csv"
    ),
    delim = ",", col_names = TRUE
    )
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
        quality_met = purrr::pmap(
          list(selected_mtx, mtx_gene_list, mtx_genefull_list),
          ~ quality_metrics(
            mtx = ..1,
            mitogenes = genes[["mito_genes"]],
            ribogenes = genes[["ribo_genes"]],
            pcgenes = genes[["protein_genes"]],
            mtx_gene = ..2,
            mtx_genefull = ..3
          )
        )
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate("Unmapped_reads_%" = sum(dplyr::c_across(tidyselect::starts_with("%_of_reads_unmapped_")), na.rm = TRUE)) %>%
      dplyr::ungroup() |>
      dplyr::select(-selected_mtx) |>
      dplyr::relocate(rna_count, quality_met, .after = library_prep)

    ### Save quality control metrics
    results <- results |>
      tidyr::unnest(quality_met)

    vroom::vroom_write(results, base::paste0(
      here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics/"),
      study_name,
      "_quality_metrics.csv"
    ),
    delim = ",", col_names = TRUE
    )
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
      dplyr::mutate(
        barcode = ic_id,
        rna_count = rna_count
      )

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
      dplyr::rename(
        nCounts = nUMIs,
        lnCounts = lnUMIs
      ) |>
      dplyr::relocate(rna_count,
        barcode,
        nCounts,
        lnCounts,
        lnFeatures,
        mitochondrial_fraction,
        ribosomal_fraction,
        coding_fraction,
        contrast_fraction,
        complexity,
        .after = library_prep
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate("Unmapped_reads_%" = sum(dplyr::c_across(tidyselect::starts_with("%_of_reads_unmapped_")), na.rm = TRUE)) %>%
      dplyr::ungroup()

    ## Save CSV
    vroom::vroom_write(
      results,
      base::paste0(
        here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics/"),
        study_name,
        "_quality_metrics.csv"
      ),
      delim = ",",
      col_names = TRUE
    )
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
        trimws(whitespace = '["; ]+', which = "left") |>
        str_extract('(?<=")[^"]+(?=")')
      tibble("{att}" := extracted)
    })
  })

  # Combine results with the original data (removing the attribute column)
  messy <- bind_cols(messy, results) |> select(-attribute)

  return(messy)
}



# Thresholds --------------------------------------------------------------
#' Calculate which barcodes that fail qc with specific thresholds
#'
#' @param qc_met A data frame with quality metrics
#' @param thresholds A data frame with thresholds
#' @param qc_droplet A vector of column names from quality control metrics
#' @param qc_plate A vector of column names from quality control metrics
#'
#' @returns A dataframe, 0 = passed, 1 = failed
check_thresholds <- function(qc_met, thresholds = NULL,
                             qc_droplet = c("nUMIs", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction", "complexity"),
                             qc_plate = c("Uniquely_mapped_reads_%", "Unmapped_reads_%", "nCounts", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction", "complexity")) {
  if (base::unique(qc_met$platform) %in% c("droplet", "plate_barcode")) {
    cols <- qc_droplet
  } else if (base::unique(qc_met$platform) %in% c("plate")) {
    cols <- qc_plate
  }

  # Pivot thresholds to long format
  thres_long <- thresholds |>
    tidyr::pivot_longer(
      cols = contains(cols),
      names_to = c("metric", "bound"),
      names_pattern = "threshold_(.*)_(lower|upper)",
      values_to = "threshold"
    ) |>
    tidyr::pivot_wider(names_from = bound, values_from = threshold) |>
    dplyr::select(-tidyselect::starts_with("threshold"))

  # Pivot qc data to long format

  qc_long <- qc_met |>
    dplyr::select(name, ic_id, ic_id_sample, ic_id_donor, ic_id_study, platform, barcode, tidyselect::all_of(cols)) |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(cols),
      names_to = "metric",
      values_to = "value"
    )

  # Join QC and thresholds
  qc_joined <- qc_long |>
    dplyr::left_join(thres_long)

  # Apply threshold logic
  # 1 = failed
  # 0 = passed

  qc_checked <- qc_joined |>
    dplyr::mutate(
      failed = dplyr::case_when(
        !is.na(lower) & value < lower ~ 1,
        !is.na(upper) & value > upper ~ 1,
        TRUE ~ 0
      )
    )

  # Pivot wide and summarise failed samples
  qc_wide <- qc_checked |>
    dplyr::select(-c(value, lower, upper)) |>
    tidyr::pivot_wider(names_from = metric, values_from = failed, names_glue = "{metric}_thres") |>
    dplyr::mutate(
      dplyr::across(tidyselect::ends_with("_thres"), ~ ifelse(is.na(.), 1, .))
    ) |>
    dplyr::mutate(
      sum_failed_thres = base::rowSums(dplyr::pick(tidyselect::ends_with("_thres")), na.rm = TRUE)
    ) |>
    dplyr::mutate(
      dplyr::across(
        tidyselect::ends_with("_thres"),
        ~ dplyr::case_when(
          . == 1 & sum_failed_thres > 1 ~ 1,
          TRUE ~ 0
        ),
        .names = "{.col}_multipass"
      )
    ) |>
    dplyr::select(-sum_failed_thres_multipass)

  return(qc_wide)
}

# plotting functions ------------------------------------------------------
my_theme <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6, color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.text = element_text(size = 4),
      legend.title = element_text(size = 4),
      legend.key.size = unit(0.3, "cm"),
      strip.text = element_text(size=5),
      plot.caption = element_text(size = 4, hjust = 0.5),
      plot.subtitle = element_text(size = 6, hjust = 0.5),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA)
    )
}


my_theme_void <- function() {
  ggplot2::theme_void() +
    ggplot2::theme(
      strip.background = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5),
      legend.text = element_text(size = 4),
      legend.title = element_text(size = 4),
      legend.key.size = unit(0.3, "cm"),
      strip.text = element_text(size=5),
      plot.caption = element_text(size = 4, hjust = 0.5),
      plot.subtitle = element_text(size = 6, hjust = 0.5),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA)
    )
}

plot_rank <- function(rank, threshold, title, non_empty) {
  p <- rank |>
    dplyr::mutate(
      not_empty = dplyr::case_when(barcode %in% non_empty$not_empty ~ "Not empty",
        .default = "Empty"
      ),
      not_empty = factor(not_empty, levels = c("Not empty", "Empty"))
    ) |>
    ggplot2::ggplot(aes(x = rank, y = counts, color = not_empty)) +
    ggrastr::geom_point_rast(size = 0.1, shape = 20, raster.dpi = 500) +
    ggplot2::scale_color_manual(values = c("Not empty" = "black", "Empty" = "#D3D3D3")) +
    ggplot2::scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::labs(
      x = "Log Rank",
      y = "Log UMIs",
      color = "Droplet",
      linetype = "Lower threshold",
      title = title
    ) +
    ggplot2::geom_hline(aes(yintercept = threshold$automatic, linetype = "Automatic"), color = "red") +
    ggplot2::scale_linetype_manual(values = c("Automatic" = 3, "Manual" = 5)) +
    my_theme() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "inside",
      legend.position.inside = c(0.2, 0.2),
      legend.key.size = unit(0.1, "inch"),
      legend.spacing.y = unit(0, "inch"),
      legend.margin = margin(0, 0, 0, 0),
      legend.text = element_text(margin = margin(l = 0)),
      legend.title = element_blank()
    )

  if (!is.null(threshold$manual)) {
    p <- p + ggplot2::geom_hline(aes(yintercept = threshold$manual, linetype = "Manual"), color = "blue")
  }

  return(p)
}

plot_rank_line <- function(rank, title, failed) {
  p <- rank |>
    dplyr::mutate(Passed = dplyr::case_when(name %in% failed_samples ~ "Failed",
      .default = "Passed"
    )) |>
    ggplot2::ggplot(aes(x = rank, y = counts, group = name, color = Passed)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = c("Passed" = "black", "Failed" = "red")) +
    ggplot2::scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::labs(
      x = "Log Rank",
      y = "Log UMIs",
      title = title
    ) +
    my_theme() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "inside",
      legend.position.inside = c(0.2, 0.2),
      legend.key.size = unit(0.1, "inch"),
      legend.spacing.y = unit(0, "inch"),
      legend.margin = margin(0, 0, 0, 0),
      legend.text = element_text(margin = margin(l = 0)),
      legend.title = element_blank()
    )

  return(p)
}

plot_hist_qc <- function(.x, .y, subtitle) {
  ggplot2::ggplot(data = tibble::tibble(value = .x), ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 100, fill = "black") +
    ggplot2::labs(
      title = dplyr::case_when(
        .y == "nUMIs" ~ "Number of unique transcripts",
        .y == "nCounts" ~ "Number of counts",
        .y == "nFeatures" ~ "Number of detected features",
        .y == "mitochondrial_fraction" ~ "Mitochondrial fraction",
        .y == "ribosomal_fraction" ~ "Ribosomal fraction",
        .y == "coding_fraction" ~ "Protein coding fraction",
        .y == "contrast_fraction" ~ "Contrast fraction",
        .y == "complexity" ~ "Library complexity",
        TRUE ~ "value",
        subtitle = subtitle
      ),
      x = dplyr::case_when(
        .y == "nUMIs" ~ "n(UMIs)",
        .y == "nCounts" ~ "n(Counts)",
        .y == "nFeatures" ~ "n(Features)",
        .y == "mitochondrial_fraction" ~ "% of counts in mitochondrial genes",
        .y == "ribosomal_fraction" ~ "% of counts in ribosomal genes",
        .y == "coding_fraction" ~ "% of counts in protein coding genes",
        .y == "contrast_fraction" ~ "exon / exon+intron counts",
        .y == "complexity" ~ "LnFeatures/Ln(Count) or ln(UMI)",
        TRUE ~ "value"
      ),
      y = "Frequency"
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) ifelse(x == 0, x, ifelse(abs(x) >= 10000 | abs(x) < 0.0001, scales::scientific(x, digits = 1), x)),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    my_theme() +
    ggplot2::theme(aspect.ratio = 1)
}

plot_hist_qc_thres <- function(.x, .y, lower, upper, subtitle = NULL, caption = NULL, .xlim = NULL) {
  ggplot2::ggplot(data = tibble::tibble(value = .x), ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 100, fill = "black") +
    ggplot2::labs(
      subtitle = subtitle,
      caption = paste(caption, "\nBlue line: Upper threshold, Red line: Lower threshold"),
      title = dplyr::case_when(
        .y == "nUMIs" ~ "Number of unique transcripts",
        .y == "nCounts" ~ "Number of counts",
        .y == "nFeatures" ~ "Number of detected features",
        .y == "mitochondrial_fraction" ~ "Mitochondrial fraction",
        .y == "ribosomal_fraction" ~ "Ribosomal fraction",
        .y == "coding_fraction" ~ "Protein coding fraction",
        .y == "contrast_fraction" ~ "Contrast fraction",
        .y == "complexity" ~ "Library complexity",
        .y == "Uniquely_mapped_reads_%" ~ "% of uniquely mapped reads",
        .y == "Unmapped_reads_%" ~ "% of unmapped reads",
        TRUE ~ "value"
      ),
      x = dplyr::case_when(
        .y == "nUMIs" ~ "n(UMIs)",
        .y == "nCounts" ~ "n(Counts)",
        .y == "nFeatures" ~ "n(Features)",
        .y == "mitochondrial_fraction" ~ "% of counts in mitochondrial genes",
        .y == "ribosomal_fraction" ~ "% of counts in ribosomal genes",
        .y == "coding_fraction" ~ "% of counts in protein coding genes",
        .y == "contrast_fraction" ~ "exon / exon+intron counts",
        .y == "complexity" ~ "LnFeatures/Ln(Count) or ln(UMI)",
        .y == "Uniquely_mapped_reads_%" ~ "% of uniquely mapped reads",
        .y == "Unmapped_reads_%" ~ "% of unmapped reads",
        TRUE ~ "value"
      ),
      y = "Frequency"
    ) +
    ggplot2::scale_x_continuous(
      limits = .xlim,
      breaks = scales::pretty_breaks(n = 5)
    ) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    ggplot2::geom_vline(xintercept = lower, color = "red", linewidth = 0.1) +
    ggplot2::geom_vline(xintercept = upper, color = "blue", linewidth = 0.1) +
    my_theme() +
    ggplot2::theme(axis.text = element_text(size = 4, color = "black")) +
    ggplot2::theme(aspect.ratio = 1)
}

plot_bar_sample_donor <- function(df, x = x, y = y, title = NULL, subtitle = NULL) {
  df |>
    ggplot2::ggplot(aes(x = x, y = y)) +
    ggplot2::geom_bar(
      stat = "identity",
      position = position_dodge(),
      fill = "black"
    ) +
    ggplot2::labs(
      subtitle = subtitle,
      title = dplyr::case_when(
        title == "nUMIs" ~ "Number of unique transcripts",
        title == "nCounts" ~ "Number of counts",
        title == "nFeatures" ~ "Number of detected features",
        title == "mitochondrial_fraction" ~ "Mitochondrial fraction",
        title == "ribosomal_fraction" ~ "Ribosomal fraction",
        title == "coding_fraction" ~ "Protein coding fraction",
        title == "contrast_fraction" ~ "Contrast fraction",
        title == "complexity" ~ "Library complexity",
        title == "Uniquely_mapped_reads_%" ~ "% of uniquely mapped reads",
        title == "Unmapped_reads_%" ~ "% of unmapped reads",
        TRUE ~ "value"
      ),
      x = dplyr::case_when(
        title == "nUMIs" ~ "n(UMIs)",
        title == "nCounts" ~ "n(Counts)",
        title == "nFeatures" ~ "n(Features)",
        title == "mitochondrial_fraction" ~ "% of counts in mitochondrial genes",
        title == "ribosomal_fraction" ~ "% of counts in ribosomal genes",
        title == "coding_fraction" ~ "% of counts in protein coding genes",
        title == "contrast_fraction" ~ "exon / exon+intron counts",
        title == "complexity" ~ "LnFeatures/Ln(Count) or ln(UMI)",
        title == "Uniquely_mapped_reads_%" ~ "% of uniquely mapped reads",
        title == "Unmapped_reads_%" ~ "% of unmapped reads",
        TRUE ~ "value"
      ),
      y = "% of total barcodes"
    ) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    my_theme() +
    ggplot2::theme(axis.text = element_text(size = 4, color = "black")) +
    ggplot2::theme(aspect.ratio = 1)
}

# Here are the groups I used to organize the metrics:
#
#   Read Counts and Percentages:
#
#   Number of input reads
# % of uniquely mapped reads
# % of chimeric reads
# % of reads mapped to multiple loci
# % of reads mapped to too many loci
# % of reads unmapped (other)
# % of reads unmapped (too many mismatches)
# % of reads unmapped (too short)
# Number of uniquely mapped reads
# Read Lengths:
#
#   Average input read length
# Average mapped read length
# Deletion and Insertion Metrics:
#
#   Average deletion length
# Deletion rate per base (%)
# Average insertion length
# Insertion rate per base (%)
# Mismatch Metrics:
#
#   Mismatch rate per base (%)
# Gene Read Counts:
#
#   Number of reads (gene sum)
# Number of reads (gene full sum)
# Number of chimeric reads
# Splice Metrics:
#
#   Number of splices (AT/AC)
# Number of splices (Annotated)
# Number of splices (GC/AG)
# Number of splices (GT/AG)
# Number of splices (Non-canonical)
# Total number of splices

plot_hist_star <- function(.x, .y) {
  ggplot2::ggplot(data = tibble::tibble(value = .x), ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 100, fill = "black") +
    ggplot2::labs(
      title = dplyr::case_when(
        .y == "Number_of_input_reads" ~ "Number of input reads",
        .y == "Uniquely_mapped_reads_%" ~ "% of uniquely mapped reads",
        .y == "Unmapped_reads_%" ~ "% of unmapped reads",
        .y == "Average_input_read_length" ~ "Average input read length",
        .y == "Average_mapped_length" ~ "Average mapped read length",
        .y == "%_of_chimeric_reads" ~ "% of chimeric reads",
        .y == "%_of_reads_mapped_to_multiple_loci" ~ "% of reads mapped to multiple loci",
        .y == "%_of_reads_mapped_to_too_many_loci" ~ "% of reads mapped to too many loci",
        .y == "%_of_reads_unmapped_other" ~ "% of reads unmapped (other)",
        .y == "%_of_reads_unmapped_too_many_mismatches" ~ "% of reads unmapped (too many mismatches)",
        .y == "%_of_reads_unmapped_too_short" ~ "% of reads unmapped (too short)",
        .y == "Deletion_average_length" ~ "Average deletion length",
        .y == "Deletion_rate_per_base_%" ~ "Deletion rate per base (%)",
        .y == "Insertion_average_length" ~ "Average insertion length",
        .y == "Insertion_rate_per_base_%" ~ "Insertion rate per base (%)",
        .y == "Mismatch_rate_per_base__%" ~ "Mismatch rate per base (%)",
        .y == "Number_of_chimeric_reads" ~ "Number of chimeric reads",
        .y == "Number_of_reads_mapped_to_multiple_loci" ~ "Number of reads mapped to multiple loci",
        .y == "Number_of_reads_mapped_to_too_many_loci" ~ "Number of reads mapped to too many loci",
        .y == "Number_of_reads_unmapped_other" ~ "Number of reads unmapped (other)",
        .y == "Number_of_reads_unmapped_too_many_mismatches" ~ "Number of reads unmapped (too many mismatches)",
        .y == "Number_of_reads_unmapped_too_short" ~ "Number of reads unmapped (too short)",
        .y == "Number_of_splices_AT/AC" ~ "Number of splices (AT/AC)",
        .y == "Number_of_splices_Annotated_(sjdb)" ~ "Number of splices (Annotated)",
        .y == "Number_of_splices_GC/AG" ~ "Number of splices (GC/AG)",
        .y == "Number_of_splices_GT/AG" ~ "Number of splices (GT/AG)",
        .y == "Number_of_splices_Non-canonical" ~ "Number of splices (Non-canonical)",
        .y == "Number_of_splices_Total" ~ "Total number of splices",
        .y == "Uniquely_mapped_reads_number" ~ "Number of uniquely mapped reads",
        TRUE ~ "value"
      ),
      x = dplyr::case_when(
        .y == "Number_of_input_reads" ~ "Number of input reads",
        .y == "Uniquely_mapped_reads_%" ~ "% of uniquely mapped reads",
        .y == "Unmapped_reads_%" ~ "% of unmapped reads",
        .y == "Average_input_read_length" ~ "Average input read length",
        .y == "Average_mapped_length" ~ "Average mapped read length",
        .y == "%_of_chimeric_reads" ~ "% of chimeric reads",
        .y == "%_of_reads_mapped_to_multiple_loci" ~ "% of reads mapped to multiple loci",
        .y == "%_of_reads_mapped_to_too_many_loci" ~ "% of reads mapped to too many loci",
        .y == "%_of_reads_unmapped_other" ~ "% of reads unmapped (other)",
        .y == "%_of_reads_unmapped_too_many_mismatches" ~ "% of reads unmapped (too many mismatches)",
        .y == "%_of_reads_unmapped_too_short" ~ "% of reads unmapped (too short)",
        .y == "Deletion_average_length" ~ "Average deletion length",
        .y == "Deletion_rate_per_base_%" ~ "Deletion rate per base (%)",
        .y == "Insertion_average_length" ~ "Average insertion length",
        .y == "Insertion_rate_per_base_%" ~ "Insertion rate per base (%)",
        .y == "Mismatch_rate_per_base__%" ~ "Mismatch rate per base (%)",
        .y == "Number_of_chimeric_reads" ~ "Number of chimeric reads",
        .y == "Number_of_reads_mapped_to_multiple_loci" ~ "Number of reads mapped to multiple loci",
        .y == "Number_of_reads_mapped_to_too_many_loci" ~ "Number of reads mapped to too many loci",
        .y == "Number_of_reads_unmapped_other" ~ "Number of reads unmapped (other)",
        .y == "Number_of_reads_unmapped_too_many_mismatches" ~ "Number of reads unmapped (too many mismatches)",
        .y == "Number_of_reads_unmapped_too_short" ~ "Number of reads unmapped (too short)",
        .y == "Number_of_splices_AT/AC" ~ "Number of splices (AT/AC)",
        .y == "Number_of_splices_Annotated_(sjdb)" ~ "Number of splices (Annotated)",
        .y == "Number_of_splices_GC/AG" ~ "Number of splices (GC/AG)",
        .y == "Number_of_splices_GT/AG" ~ "Number of splices (GT/AG)",
        .y == "Number_of_splices_Non-canonical" ~ "Number of splices (Non-canonical)",
        .y == "Number_of_splices_Total" ~ "Total number of splices",
        .y == "Uniquely_mapped_reads_number" ~ "Number of uniquely mapped reads"
      ),
      y = "Frequency"
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) ifelse(x == 0, x, ifelse(abs(x) >= 10000 | abs(x) < 0.0001, scales::scientific(x, digits = 1), x)),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    my_theme() +
    ggplot2::theme(aspect.ratio = 1)
}


# Specific plots ----------------------------------------------------------
visualize_qc_hist_with_thresholds <- function(qc_met, thresholds, checked, subtitle = NULL,
                                              qc_droplet = c("nUMIs", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction", "complexity"),
                                              qc_plate = c("Uniquely_mapped_reads_%", "Unmapped_reads_%", "nCounts", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction")) {
  if (base::unique(qc_met$platform) %in% c("droplet", "plate_barcode")) {
    cols <- qc_droplet
  } else if (base::unique(qc_met$platform) %in% c("plate")) {
    cols <- qc_plate
  }

  # Number of barcodes that fail quality control
  df_thres <- checked |>
    dplyr::select(tidyselect::all_of(base::paste0(cols, "_thres"))) |>
    dplyr::summarise_all(sum) |>
    dplyr::rename_all(~ stringr::str_replace(., "_thres", ""))

  # Number of barcodes that fail more than one quality control
  df_multi <- checked |>
    dplyr::select(dplyr::ends_with("thres_multipass")) |>
    dplyr::summarise_all(sum) |>
    dplyr::rename_all(~ stringr::str_replace(., "_thres_multipass", ""))

  # Total number of barcodes
  n_total <- checked |> base::nrow()

  qc_met |>
    dplyr::select(tidyselect::all_of(cols)) |>
    purrr::iwalk(~ {
      lower_col <- paste0("threshold_", .y, "_lower")
      upper_col <- paste0("threshold_", .y, "_upper")
      lower <- thresholds[[lower_col]][1]
      upper <- thresholds[[upper_col]][1]

      # Check if lower or upper is NA or NULL and replace accordingly
      if (base::is.na(lower) || base::is.null(lower)) {
        lower <- NULL
      }

      if (base::is.na(upper) || base::is.null(upper)) {
        upper <- NULL
      }

      # How many barcodes fail threshold
      col <- base::paste0(.y)
      fail <- df_thres[[col]][1]
      fail_multi <- df_multi[[col]][1]

      # Create caption
      caption <- base::paste0("Total: ", n_total, ", Failed: ", fail, ", Failed multi: ", round((fail_multi / fail) * 100, 2), "%")

      # Plot
      base::print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle, caption = caption))
    })
}

visualize_barcodes_failed_qc <- function(checked, subtitle = NULL,
                                         qc_droplet = c("nUMIs", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction", "complexity"),
                                         qc_plate = c("Uniquely_mapped_reads_%", "Unmapped_reads_%", "nCounts", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction")) {
  
  if (base::unique(checked$platform) %in% c("droplet", "plate_barcode")) {
    cols <- qc_droplet

    # Number of barcodes that fail quality control
    df_thres <- checked |>
      dplyr::select(ic_id_sample, tidyselect::all_of(base::paste0(cols, "_thres"))) |>
      dplyr::group_by(ic_id_sample) |> 
      dplyr::summarise_all(sum) |>
      dplyr::rename_all(~ stringr::str_replace(., "_thres", ""))

    # Total number of barcodes per sample
    n_total <- checked|>
      dplyr::group_by(ic_id_sample) |>
      dplyr::tally() |>
      dplyr::full_join(df_thres, by = "ic_id_sample") |>
      dplyr::mutate(
        dplyr::across(
          c(-ic_id_sample, -n),
          .fns = function(qc) {
            round((qc / n) * 100, 2)
          }
        ),
        ic_id_sample = as.character(ic_id_sample)
      ) |>
      dplyr::select(-n) |>
      tidyr::pivot_longer(-ic_id_sample, names_to = "qc", values_to = "value")

    # Identify rows with missing values
    missing_values <- n_total |> dplyr::filter(is.na(value))
    print("Study:")
    print(subtitle)
    print("Rows with missing values:")
    print(missing_values)

    # Identify rows with values outside the scale range (assuming scale range is 0-100)
    outside_scale_range <- n_total |> dplyr::filter(value < 0 | value > 100)
    print("Rows with values outside the scale range:")
    print(outside_scale_range)

    plot <- n_total |>
      ggplot2::ggplot(aes(y = ic_id_sample, x = value)) +
      ggplot2::geom_bar(
        stat = "identity",
        position = position_dodge(),
        fill = "black"
      ) +
      ggplot2::geom_text(aes(label = value), size = 1.5, hjust = -0.2) +
      ggplot2::labs(
        subtitle = subtitle,
        x = "% of total barcodes",
        y = "Sample"
      ) +
      ggplot2::scale_x_continuous(limits = c(0, 100), breaks = scales::pretty_breaks(n = 5)) +
      ggplot2::facet_wrap(dplyr::case_when(
        qc == "nUMIs" ~ "Number of unique transcripts",
        qc == "nCounts" ~ "Number of counts",
        qc == "nFeatures" ~ "Number of detected features",
        qc == "mitochondrial_fraction" ~ "Mitochondrial fraction",
        qc == "ribosomal_fraction" ~ "Ribosomal fraction",
        qc == "coding_fraction" ~ "Protein coding fraction",
        qc == "contrast_fraction" ~ "Contrast fraction",
        qc == "complexity" ~ "Library complexity",
        qc == "Uniquely_mapped_reads_%" ~ "% of uniquely mapped reads",
        qc == "Unmapped_reads_%" ~ "% of unmapped reads",
        TRUE ~ qc
      ) ~ ., nrow = 1) +
      my_theme()

    print(plot)
  } else if (base::unique(checked$platform) %in% c("plate")) {
    cols <- qc_plate

    df_thres <- checked |>
      dplyr::select(ic_id_donor, tidyselect::all_of(base::paste0(cols, "_thres"))) |>
      dplyr::group_by(ic_id_donor) |> 
      dplyr::summarise_all(sum) |>
      dplyr::rename_all(~ stringr::str_replace(., "_thres", ""))

    # Total number of barcodes
    n_total <- checked |>
      dplyr::group_by(ic_id_donor) |>
      dplyr::tally() |>
      dplyr::full_join(df_thres, by = "ic_id_donor") |>
      dplyr::mutate(
        dplyr::across(
          c(-ic_id_donor, -n),
          .fns = function(qc) {
            round((qc / n) * 100, 2)
          }
        ),
        ic_id_donor = as.character(ic_id_donor)
      ) |>
      dplyr::select(-n) |>
      tidyr::pivot_longer(-ic_id_donor,
        names_to = "qc",
        values_to = "value"
      )

    # Identify rows with missing values
    missing_values <- n_total |> dplyr::filter(is.na(value))
    print("Study:")
    print(subtitle)
    print("Rows with missing values:")
    print(missing_values)

    # Identify rows with values outside the scale range (assuming scale range is 0-100)
    outside_scale_range <- n_total |> dplyr::filter(value < 0 | value > 100)
    print("Rows with values outside the scale range:")
    print(outside_scale_range)

    plot <- n_total |>
      ggplot2::ggplot(aes(y = ic_id_donor, x = value)) +
      ggplot2::geom_bar(
        stat = "identity",
        position = position_dodge(),
        fill = "black"
      ) +
      ggplot2::geom_text(aes(label = value), size = 1.5, hjust = -0.2) +
      ggplot2::labs(
        subtitle = subtitle,
        x = "% of total barcodes",
        y = "Donor"
      ) +
      ggplot2::scale_x_continuous(limits = c(0, 100), breaks = scales::pretty_breaks(n = 5)) +
      ggplot2::facet_wrap(dplyr::case_when(
        qc == "nUMIs" ~ "Number of unique transcripts",
        qc == "nCounts" ~ "Number of counts",
        qc == "nFeatures" ~ "Number of detected features",
        qc == "mitochondrial_fraction" ~ "Mitochondrial fraction",
        qc == "ribosomal_fraction" ~ "Ribosomal fraction",
        qc == "coding_fraction" ~ "Protein coding fraction",
        qc == "contrast_fraction" ~ "Contrast fraction",
        qc == "complexity" ~ "Library complexity",
        qc == "Uniquely_mapped_reads_%" ~ "% of uniquely mapped reads",
        qc == "Unmapped_reads_%" ~ "% of unmapped reads",
        TRUE ~ qc
      ) ~ ., nrow = 1) +
      my_theme()

    print(plot)
  }
}

visualize_qc_hist_threshold_per_sample <- function(qc_met_split, subtitle = NULL, 
                                                   qc_droplet = c("nUMIs", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction", "complexity"),
                                                   qc_plate = c("Uniquely_mapped_reads_%", "Unmapped_reads_%", "nCounts", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction")) {
  
  if (base::unique(qc_met_split$platform) %in% c("droplet", "plate_barcode")) {
    metric_cols <- qc_droplet
    
    qc_met_split |>
      dplyr::select(dplyr::all_of(metric_cols)) |>
      purrr::imap(~ {
        lower_col <- paste0("threshold_", .y, "_lower")
        upper_col <- paste0("threshold_", .y, "_upper")
        
        lower <- qc_met_split[[lower_col]][1]
        upper <- qc_met_split[[upper_col]][1]
        
        if (is.na(lower) || is.null(lower)) lower <- NULL
        if (is.na(upper) || is.null(upper)) upper <- NULL
        
        print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle))
      })
  } else if (base::unique(qc_met_split$platform) %in% c("plate")) {
    metric_cols <- qc_plate
    
    qc_met_split |>
      dplyr::select(dplyr::all_of(metric_cols)) |>
      purrr::imap(~ {
        lower_col <- paste0("threshold_", .y, "_lower")
        upper_col <- paste0("threshold_", .y, "_upper")
        
        lower <- qc_met_split[[lower_col]][1]
        upper <- qc_met_split[[upper_col]][1]
        
        if (is.na(lower) || is.null(lower)) lower <- NULL
        if (is.na(upper) || is.null(upper)) upper <- NULL
        
        print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle)) 
        
      }) 
  }
}


# Converting --------------------------------------------------------------
parse_to_hours <- function(time_str) {
  parts <- strsplit(time_str, " ")[[1]]
  days <- as.numeric(gsub(" days", "", parts[1]))
  time_parts <- strsplit(parts[3], ":")[[1]]
  hours <- as.numeric(time_parts[1])
  minutes <- as.numeric(time_parts[2])
  seconds <- as.numeric(time_parts[3])
  return((days * 24) + (hours) + (minutes / 60) + seconds / 3600)
}


# barcode file processing -------------------------------------------------

#' Add prefix to barcode or replace barcode with ic_id
#'
#' @param exists bool; TRUE/ FALSE does path to barcode file exsist with file.exsist()
#' @param path Path to barcode file
#' @param ic_id_sample Name to add as prefix or replace barcode with
#' @param platform A string (plate, plate_barcode or droplet) will determine how ic_id_sample is handled
#' @param save_path Path to save the new barcode file
#'
#' @returns
#' @export
#'
#' @examples
process_barcode_files <- function(exists, path, ic_id_sample, platform, save_path){
  
  # If file does not exist, stop
  if (!exists){
    stop("path does not exist!")
  } 
  
  result <-  tryCatch({
    # Load barcode file 
    barcodes <- vroom::vroom_lines(path)
    
    # Replace or prefix barcode
    if (platform %in% c("droplet", "plate_barcode")){
      prefixed_lines <- base::paste0(ic_id_sample, "_", barcodes)
    } else if (platform %in% c("plate")) {
      prefixed_lines <- ic_id_sample 
    } else {
      stop("Platform does not match any of the requirements!")
    }
    
    # Save barcode file 
    vroom::vroom_write_lines(x = prefixed_lines, file = save_path)
    
    # Check that the new saved file actually exists
    if (!file.exists(save_path)) {
      stop(paste0("Did not successfully write to: ", save_path))
    }
  },
  error = function(e) {
    message("Error processing ", path, ": ", e$message)
  })
}


# GO-term -----------------------------------------------------------------
go_term_wald <- function(test_genes, bg_genes, bg_cell){
  
  # add entrez IDs to genes
  
  genes <- BiocGenerics::Map(function(df){
    if (length(rownames(df) > 0)){
      df$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                        keys=df$gene,
                                        column="ENTREZID",
                                        keytype="SYMBOL",
                                        multiVals="first")
      
      output <- df
    } else {
      output <- "NA"
    }
    
    return(output)
  },
  df = test_genes)
  
  # Add entrez IDs to background genes
  
  bg_genes <- BiocGenerics::Map(function(df){
    df$entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                      keys=df$gene,
                                      column="ENTREZID",
                                      keytype="SYMBOL",
                                      multiVals="first")
    
    return(df)
  },
  df = bg_genes)
  
  # Create universe - this should be any gene that could have been positive
  # And I guess that should be all genes that were tested
  
  my_universe <- BiocGenerics::Map(function(df){
    output <- df %>%
      as.data.frame() %>%
      filter(!is.na(entrez)) %>%
      pull(entrez) %>%
      unlist() %>%
      unname() %>%
      unique()
    
    return(output)
  },
  df = bg_genes)
  
  gene_list <- genes
  # Pathway analysis with clusterprofiler + jaccard index
  
  # Jaccard
  # k is the overlap between your genes-of-interest and the geneset
  # n is the number of all unique genes-of-interest
  
  # BgRatio=M/N
  
  # M is the number of genes within each geneset
  # N is the number of all unique genes across all genesets (universe)
  output <- BiocGenerics::Map(function(df, universe){
    if (is.character(df) == TRUE) {
      output <- "no genes"
    } else {
      output <- df %>%
        dplyr::select(entrez) %>%
        purrr::as_vector() %>%
        unname() %>%
        unique() %>%
        clusterProfiler::enrichGO(OrgDb = "org.Hs.eg.db",
                                  pAdjustMethod = "fdr",
                                  ont = "BP",
                                  pvalueCutoff  = 0.2,
                                  qvalueCutoff  = 0.2,
                                  readable      = TRUE,
                                  universe = universe
        ) %>%
        dplyr::mutate(is_sig = case_when(p.adjust <= 0.01 ~ "*",
                                         p.adjust > 0.01 ~ "ns"))
    }
    
    return(output)
    
  },
  df = gene_list,
  universe = my_universe)
  
  
  return(output)
}

