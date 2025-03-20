my_theme <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(strip.background = element_blank(),
                   plot.title = element_text(size = 8),
                   axis.title = element_text(size = 7),
                   axis.text = element_text(size = 6))
}

#' Create Directories if They Do Not Exist
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
                  name = i) |>
    dplyr::relocate(name)
  return(df)
}

read_mtx <- function(path, feature, workers =  parallelly::availableCores()){
  
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
  
  # Read files in parallel
  future::plan(multisession, workers = workers)  # Enable parallel processing
  
  mtx_list <- furrr::future_pmap(file_group, function(mtx, barcodes, features) {
    Seurat::ReadMtx(mtx = mtx, cells = barcodes, features = features, feature.column = 1)
  },
  .options = furrr::furrr_options(seed = TRUE))|> 
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

prefix_colnames <- function(mtx_list, workers = parallelly::availableCores()) {
  # Enable parallel processing
  future::plan(multisession, workers = workers)
  
  furrr::future_imap(mtx_list, ~ {
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
  }, .options = furrr::furrr_options(seed = TRUE)) 
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

rank_plot <- function(rank, mtx_filter){
  rank$ranks |> 
    tibble::rownames_to_column("barcode") |> 
    dplyr::mutate(not_empty = dplyr::case_when(barcode %in% colnames(mtx_filter) ~ "Not empty",
                                               .default = base::as.character("Empty")),
                  not_empty = base::factor(not_empty, levels = c("Not empty", "Empty"))) |> 
    ggplot2::ggplot(aes(x = rank, y = counts, color = not_empty)) +
    ggplot2::geom_point(shape = 20) +
    ggplot2::scale_color_manual(values = c("Not empty" = "#000000", "Empty" = "#D3D3D3")) +
    ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    ggplot2::labs(x = "Log Rank",
                  y = "Log UMIs") +
    ggplot2::geom_hline(yintercept = rank$lower.threshold, color = "red", linetype = "dashed") +
    my_theme()
}

remove_empty_droplets <- function(mtx, fdr = 0.001){
  # Get barcode ranks and threshold
  rank<- rank_barcodes(mtx)
  
  # Find empty droplets
  empty <- DropletUtils::emptyDrops(mtx, 
                                    lower = rank$lower.threshold, 
                                    BPPARAM = BiocParallel::MulticoreParam())
  
  # Remove empty droplets
  mtx_filter <-
    mtx[, BiocGenerics::colnames(mtx) %in%
          BiocGenerics::rownames(empty)[base::which(empty$FDR <= fdr)]]
  
  # Rank plot
  plot <- rank_plot(rank = rank, mtx_filter = mtx_filter)
  
  output <- base::list(rank = rank$ranks,
                       lower_threshold = rank$lower.threshold,
                       not_empty = BiocGenerics::colnames(mtx_filter),
                       rank_plot= plot)
}

quality_metrics <- function(mtx, barcodes, mitogenes, ribogenes, pcgenes, mtx_gene, mtx_genefull){
  # keep only non empty droplets
  none_mtx <- mtx[, BiocGenerics::colnames(mtx) %in% barcodes]
  
  ## Contrast 
  # Match barcodes between the input counts and the mtx_gene matrix
  mtx_gene_sub <- mtx_gene[, BiocGenerics::colnames(mtx_gene) %in% BiocGenerics::colnames(none_mtx)]
  mtx_gene_sub <- mtx_gene_sub[, base::match(BiocGenerics::colnames(none_mtx), BiocGenerics::colnames(mtx_gene_sub))]
  # Match barcodes between the input counts and the mtx_gene matrix
  mtx_genefull_sub <- mtx_genefull[, BiocGenerics::colnames(mtx_genefull) %in% BiocGenerics::colnames(none_mtx)]
  mtx_genefull_sub <- mtx_genefull_sub[, base::match(BiocGenerics::colnames(none_mtx), BiocGenerics::colnames(mtx_genefull_sub))]

  # Internalize matrix
  metrics <- BiocGenerics::as.data.frame(base::matrix(ncol=8, 
                                                        nrow=ncol(none_mtx)))
  # Set colnames
  BiocGenerics::colnames(metrics) <- c("barcode",
                                       "logUMIs",
                                       "logFeatures",
                                       "mitochondrial_fraction",
                                       "ribosomal_fraction",
                                       "coding_fraction", 
                                       "contrast_fraction", 
                                       "complexity")
  
  metrics[,1] <- BiocGenerics::colnames(none_mtx)
  metrics[,2] <- Matrix::colSums(none_mtx)
  metrics[,3] <- Matrix::colSums(none_mtx > 0)
  metrics[,4] <- Matrix::colSums(none_mtx[ BiocGenerics::rownames(none_mtx) %in% mitogenes,]) / metrics[,2]
  metrics[,5] <- Matrix::colSums(none_mtx[ BiocGenerics::rownames(none_mtx) %in% ribogenes,]) / metrics[,2]
  metrics[,6] <- Matrix::colSums(none_mtx[ BiocGenerics::rownames(none_mtx) %in% pcgenes,]) / metrics[,2]
  metrics[,7] <- Matrix::colSums(mtx_gene_sub) / Matrix::colSums(mtx_genefull_sub)
  metrics[,8] <- metrics[,3] / metrics[,2]
  metrics[,2] <- base::log(metrics[,2])
  metrics[,3] <- base::log(metrics[,3])
  
  return(metrics)
}

return(mtx_list)
process_study_samples <- function(path, study_metadata, genes, study_name, droplet_based, plate_based) {

  
  # Extract unique library prep
  library_prep <- base::unique(study_metadata$library_prep)
  if (length(library_prep) > 1) stop("Multiple library preps found in study. Please split them first.")
  
  # Check if it's droplet or plate-based
  if (library_prep %in% droplet_based) {
    base::message("Method is droplet based")
    
    # Load all mtx files
    mtx_gene_list <- read_mtx(path, "Gene")
    mtx_gene_list <- sample_to_ic_id(mtx_gene_list, study_metadata) # Add ic_id names to mtx list
    mtx_gene_list <- prefix_colnames(mtx_gene_list) # Add ic_id prefix to barcodes
    
    mtx_genefull_list <- read_mtx(path, "GeneFull")
    mtx_genefull_list <- sample_to_ic_id(mtx_genefull_list, study_metadata)
    mtx_genefull_list <- prefix_colnames(mtx_genefull_list)
    
    # Load GeneFull_Ex50pAS only if any sample is nuclei
    mtx_genefull_50_list <- NULL
    if ("nuclei" %in% study_metadata$cell_nuclei) {
      base::message("Nuclei samples, loading GeneFull_Ex50pAS")
      
      mtx_genefull_50_list <- read_mtx(path, "GeneFull_Ex50pAS")
      mtx_genefull_50_list <- sample_to_ic_id(mtx_genefull_50_list, study_metadata)
      mtx_genefull_50_list <- prefix_colnames(mtx_genefull_50_list)
    }
    
    future::plan(multisession, workers = parallelly::availableCores())  # Enable parallel processing
    
    # Process each sample individually
    results <- study_metadata %>%
      dplyr::mutate(
        # Main mtx file to use, if it is cell use gene count, if nuclei use genefull_ex50pAS counts
        selected_mtx = purrr::map2(cell_nuclei, ic_id, ~ {
          if (.x == "cell") {
            base::message(paste("Selected mtx_gene for ic_id:", .y))
            mtx_gene_list[[.y]]
          } else {
            base::message(paste("Selected mtx_genefull_50 for ic_id:", .y))
            mtx_genefull_50_list[[.y]]
          }
        })
        # Identify empty droplets
        barcodes = purrr::map(selected_mtx, remove_empty_droplets),
        # Quality control per sample
        quality_met = furrr::future_pmap(list(selected_mtx, barcodes, mtx_gene_list, mtx_genefull_list),
                                  ~ quality_metrics(
                                    mtx = ..1,
                                    barcodes = ..2[["not_empty"]],
                                    mitogenes = genes[["mito_genes"]],
                                    ribogenes = genes[["ribo_genes"]],
                                    pcgenes = genes[["protein_genes"]],
                                    mtx_gene = ..3,
                                    mtx_genefull = ..4
                                  ))) 

    
  } else if (library_prep %in% plate_based) {
    base::message("Method is plate based")
    # Load all samples into a single matrix
    mtx_gene <- read_mtx(path, "Gene")
    mtx_genefull <- read_mtx(path, "GeneFull")
    
    mtx_genefull_50 <- NULL
    if ("nuclei" %in% study_metadata$cell_nuclei) {
      base::message("Nuclei samples, loading GeneFull_Ex50pAS")
      mtx_genefull_50 <- read_mtx(path, "GeneFull_50")
    }
    
    # Use GeneFull_50 for QC if nuclei, otherwise Gene
    selected_mtx <- if ("nuclei" %in% study_metadata$cell_nuclei) mtx_genefull_50 else mtx_gene
    
    # No barcode filtering for plate-based
    results <- list(study_name = quality_metrics(
      mtx = selected_mtx,
      barcodes = colnames(selected_mtx),
      mitogenes = genes[["mito_genes"]],
      ribogenes = genes[["ribo_genes"]],
      pcgenes = genes[["protein_genes"]],
      mtx_gene = mtx_gene,
      mtx_genefull = mtx_genefull
    ))
  } else {
    stop("Unknown library preparation type")
  }
  
  return(results)
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
  
