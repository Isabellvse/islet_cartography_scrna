# Description -------------------------------------------------------------
# From quality control plots, add threshold data for each study and specific donors / samples
# We will also add thresholds and annotate whether a sample or donor is excluded

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)
create_directories(here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics_updated"))

# Load --------------------------------------------------------------------
## Quality metrics ----
paths <- base::list.files(
  path = here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics"),
  pattern = "quality_metrics.csv",
  full.names = TRUE,
  recursive = TRUE
)

base::names(paths) <- purrr::map_chr(paths, ~ stringr::str_extract(.x, pattern = "(?<=quality_metrics/)[^/]+(?=_quality_metrics\\.csv)"))

quality_met <- purrr::map(paths, ~ vroom::vroom(.x, delim = ",", col_names = TRUE)) %>%
  purrr::modify_depth(1, ~ dplyr::rowwise(.) %>%
    # Calculate % of unmapped reads
    dplyr::mutate(
      "Unmapped_reads_%" = sum(dplyr::c_across(tidyselect::starts_with("%_of_reads_unmapped_")), na.rm = TRUE),
      # Make ic_ids character
      dplyr::across(.cols = tidyselect::starts_with("ic_id_"), .fns = as.character)
    ) %>%
    dplyr::ungroup() %>%
    # Removed failed samples from barcode
    dplyr::filter(!ic_id %in% c("ic_25_11_11", "ic_25_7_7"), !donor == "excluded"))

## Thresholds
# Split into a list
thresholds <- vroom::vroom(here::here("islet_cartography_scrna/data/quality_control/first_pass/first_pass_threshold.csv"),
  delim = ";",
  col_names = TRUE
) |>
  (function(df) base::split(df, base::factor(df$name)))()

# Thresholds not split into a list
thresholds_df <- vroom::vroom(here::here("islet_cartography_scrna/data/quality_control/first_pass/first_pass_threshold.csv"),
  delim = ";",
  col_names = TRUE
)

# Preprocess --------------------------------------------------------------

# Divide Kang into cell and nuclei
quality_met[["Kang_cell"]] <- quality_met[["Kang"]] |>
  dplyr::filter(cell_nuclei == "cell") |>
  dplyr::mutate(name = paste0(name, "_cell"))
quality_met[["Kang_nuclei"]] <- quality_met[["Kang"]] |>
  dplyr::filter(cell_nuclei == "nuclei") |>
  dplyr::mutate(name = paste0(name, "_nuclei"))
quality_met[["Kang"]] <- NULL


# Add column which states is a donor / sample is excluded -----------------
quality_met <- quality_met |>
  purrr::modify_depth(1, ~ dplyr::mutate(.,
    "excluded" = dplyr::case_when(
      ic_id %in% c(
        "ic_25_11_11", # Poor rank plot
        "ic_25_7_7", # Poor rank plot
        "ic_11_12_359" # Only one cell
      ) ~ "excluded",
      donor == "excluded" ~ "excluded",
      .default = "included"
    )
  ))

# Add donor specific thresholds -------------------------------------------
sample_thresholds <- data.frame(
  "name" = c(rep("Gurp", 2), rep("HPAP_10x", 4), rep("Wang_Sander", 4)),
  "ic_id_sample" = c("5", "6", "9", "30", "13", "14", "4", "3", "16", "19"),
  "threshold_nUMIs_lower" = c(10000, 10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000),
  "threshold_nUMIs_upper" = c(90000, 50000, 50000, 15000, 50000, 50000, 10000, 10000, 10000, 10000),
  "threshold_nFeatures_lower" = c(3000, 3000, 1000, 1000, 1000, 2000, 500, 500, 500, 500),
  "threshold_nFeatures_upper" = c(9000, 8000, 7000, 7000, 7000, 9000, 10000, 10000, 10000, 10000),
  "threshold_complexity_lower" = c(0.3, 0.1, 0.4, 0.4, 0.3, 0.3, 0.4, 0.3, 0.3, 0.4),
  "threshold_complexity_upper" = c(0.8, 0.55, 1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)
)

# Update thresholds -------------------------------------------------------
thresholds_updated <- purrr::imap(quality_met, function(qc, name) {
  # Add thresholds to each donor / sample
  thresholds_original <- qc |>
    dplyr::select(name, ic_id, ic_id_sample, ic_id_donor) |>
    # Make dataframe distinct to remove repeated rows
    dplyr::distinct() |>
    dplyr::left_join(
      y = thresholds[[name]],
      by = c("name"),
      relationship = "many-to-one"
    )
  # Update rows that have new thresholds
  thresholds_updated <- thresholds_original |>
    dplyr::rows_update(y = sample_thresholds, by = c("name", "ic_id_sample"), unmatched = "ignore")

  # check that rows you wanted updated actually are
  # return all rows in x without a match in y
  not_match <- dplyr::anti_join(
    x = thresholds_updated,
    y = thresholds_original,
    by = names(thresholds_updated)[names(thresholds_updated) %in% names(thresholds_original)]
  )

  if (nrow(not_match > 0)) {
    print(not_match)
  }
  return(thresholds_updated)
})

# Add sample / donor specific thresholds to overall thresholds ------------
# 0 = passed
# 1 = failed

checked <- purrr::imap(quality_met, function(qc, name) {
  checked <- check_thresholds(qc_met = qc, thresholds = thresholds_updated[[name]])
  return(checked)
})


# Combine thresholds, passing and quality control together ----------------
# Updated thresholds
quality_met_updated <- imap(quality_met, function(df, name) {
  df |>
    dplyr::left_join(
      y = thresholds_updated[[name]],
      relationship = "many-to-one"
    ) |>
    dplyr::full_join(y = checked[[name]])
})


# Plot --------------------------------------------------------------------
## Barcodes removed per sample / donor ----
pdf(
  file = here::here(
    "islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_per_donor_sample_bar_updated.pdf"
  ),
  height = 3,
  width = 16 # Adjust width to accommodate multiple facets in one row
)
purrr::iwalk(checked, function(df, subtitle) {
  subtitle <- subtitle
  visualize_barcodes_failed_qc(checked = df, subtitle = subtitle)
})
dev.off()


## qc metrics per sample / donor ----
pdf(
  file = here::here(
    "islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_per_donor_sample_updated.pdf"
  ),
  height = 2,
  width = 2
)

# Split by donor /sample
quality_met_updated_split <- quality_met_updated |>
  purrr::modify_depth(1, ~ if (base::unique(.x$platform) %in% c("droplet", "plate_barcode")) {
    base::split(.x, .x$ic_id_sample)
  } else if (base::unique(.x$platform) %in% c("plate")) {
    base::split(.x, .x$ic_id_donor)
  }) |>
  purrr::list_flatten()

# plot qc in histogram for each donor / sample
purrr::iwalk(quality_met_updated_split, function(df, name) {
  subtitle <- name
  print(subtitle)
  visualize_qc_hist_threshold_per_sample(qc_met_split = df, subtitle = subtitle)
})

dev.off()

## Compare changed thresholds ----
# Old thresholds
quality_met_old <- imap(quality_met, function(df, name) {
  df |>
    dplyr::left_join(
      y = thresholds[[name]],
      relationship = "many-to-one"
    )
})

# samples that were modified
to_keep <- paste0(sample_thresholds$name, "_", sample_thresholds$ic_id_sample)

# Data frame for old thresholds
quality_met_old_split <- quality_met_old |>
  purrr::modify_depth(1, ~ if (base::unique(.x$platform) %in% c("droplet", "plate_barcode")) {
    base::split(.x, .x$ic_id_sample)
  } else if (base::unique(.x$platform) %in% c("plate")) {
    base::split(.x, .x$ic_id_donor)
  }) |>
  purrr::list_flatten() |>
  (function(x) x[to_keep])() |>
  purrr::set_names(function(x) paste0(x, "_old"))

# Data frame for updated thresholds
quality_met_updated_split <- quality_met_updated |>
  purrr::modify_depth(1, ~ if (base::unique(.x$platform) %in% c("droplet", "plate_barcode")) {
    base::split(.x, .x$ic_id_sample)
  } else if (base::unique(.x$platform) %in% c("plate")) {
    base::split(.x, .x$ic_id_donor)
  }) |>
  purrr::list_flatten() |>
  (function(x) x[to_keep])() |>
  purrr::set_names(function(x) paste0(x, "_updated"))

# Find columns that a different between old thresholds and updated thresholds
qc_diff_columns <- purrr::map2(
  quality_met_old_split,
  quality_met_updated_split,
  ~ {
    # Identify which columns are different
    diff_columns <- purrr::keep(names(.x), function(colname) {
      !base::identical(.x[[colname]], .y[[colname]]) # Check if column values are different
    })

    # Return a list with only the differing columns
    colnames(.x[diff_columns]) |>
      stringr::str_remove("threshold_") |>
      stringr::str_remove("_lower") |>
      stringr::str_remove("_upper") |>
      unique()
  }
)

# Flatten the differing columns and rename
qc_diff_columns_flat <- qc_diff_columns |>
  purrr::list_flatten()

# Get columns of qc metrics
qc_colnames <- c(qc_diff_columns_flat, qc_diff_columns_flat |> purrr::set_names(function(x) stringr::str_replace(x, "_old", "_updated")))

# Plot histograms with old thresholds
plots_old <- purrr::imap(quality_met_old_split, function(df, name) {
  subtitle <- name
  print(subtitle)
  plot <- visualize_qc_hist_threshold_per_sample(qc_met_split = df, subtitle = subtitle, qc_droplet = qc_colnames[[name]], qc_plate = qc_colnames[[name]])
  return(patchwork::wrap_plots(plot, nrow = 1))
})

# Plot histograms with updated thresholds
plots_updated <- purrr::imap(quality_met_updated_split, function(df, name) {
  subtitle <- name
  print(subtitle)
  plot <- visualize_qc_hist_threshold_per_sample(qc_met_split = df, subtitle = subtitle, qc_droplet = qc_colnames[[name]], qc_plate = qc_colnames[[name]])
  return(patchwork::wrap_plots(plot, nrow = 1))
})

# Save old and updated thresholds side-by-side
pdf(
  file = here::here(
    "islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_per_donor_sample_comparison.pdf"
  ),
  height = 6,
  width = 6
)

purrr::walk2(plots_old, plots_updated, function(old, updated) {
  # Extract plots either from patchwork or as single-element list
  get_plots <- function(p) {
    if (!is.null(p$patches$plots) && length(p$patches$plots) > 0) {
      return(p$patches$plots)
    } else {
      return(list(p)) # single ggplot object
    }
  }

  old_plots <- get_plots(old)
  updated_plots <- get_plots(updated)

  # Determine the max number of columns
  n_old <- length(old_plots)
  n_updated <- length(updated_plots)
  n_max <- max(n_old, n_updated, 3)

  # Pad with plot_spacer() if needed
  if (n_old < n_max) {
    old_plots <- c(old_plots, replicate(n_max - n_old, patchwork::plot_spacer(), simplify = FALSE))
  }

  if (n_updated < n_max) {
    updated_plots <- c(updated_plots, replicate(n_max - n_updated, patchwork::plot_spacer(), simplify = FALSE))
  }

  # Wrap and stack
  old_row <- patchwork::wrap_plots(old_plots, nrow = 1)
  updated_row <- patchwork::wrap_plots(updated_plots, nrow = 1)
  combined <- old_row / updated_row

  print(combined)
})

dev.off()


# Excluded barcodes -------------------------------------------------------
# barcodes with sum_failed_thres > 0 will be exlucded as they failed quality control
# and save as csv, additionally donors with less then 5 cells will be excluded
quality_met_updated <- purrr::imap(quality_met_updated, function(df, study_name) {
  # Exclude or include barcodes
  df_q <- df |>
    dplyr::mutate(excluded = dplyr::case_when(
      excluded == "excluded" ~ "excluded",
      excluded == "included" & sum_failed_thres > 0 ~ "excluded",
      excluded == "included" & sum_failed_thres == 0 ~ "included"
    ))
  
  # Calculate how many cells are excluded and included
  df_n <- df_q |> 
    dplyr::group_by(ic_id_donor, excluded) |> 
    dplyr::tally() |> 
    dplyr::ungroup() |> 
    tidyr::pivot_wider(id_cols = ic_id_donor, 
                       names_from = excluded, 
                       values_from = n, 
                       names_prefix = "n_donor_")
  
  # Join dataframes and filter based on this thresholds
  df_join <- df_q |> dplyr::left_join(y = df_n, by = "ic_id_donor")
  
  # Exclude donors which have less than 5 cells
  df_out <- df_join |> 
    dplyr::mutate(excluded = case_when(n_donor_included <= 5 ~ "excluded", 
                                               .default = as.character(excluded)))
  # Save csv files
  vroom::vroom_write(df_out, base::paste0(
    here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics_updated/"),
    study_name,
    "_quality_metrics.csv"
  ),
  delim = ",",
  col_names = TRUE
  )

  return(df_out)
})


# Number of barcodes left -------------------------------------------------
# Count number of barcodes left following filtering
excluded_included_sum <- quality_met_updated |>
  purrr::modify_depth(1, ~ if (base::unique(.x$platform) %in% c("droplet", "plate_barcode")) {
    .x |>
      dplyr::mutate(ic_id_sample = as.character(ic_id_sample)) |>
      dplyr::group_by(ic_id_sample, excluded, platform) |>
      dplyr::tally()
  } else if (base::unique(.x$platform) %in% c("plate")) {
    .x |>
      dplyr::mutate(ic_id_donor = as.character(ic_id_donor)) |>
      dplyr::group_by(ic_id_donor, excluded, platform) |>
      dplyr::tally()
  })

pdf(
  file = here::here(
    "islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_n_included_excluded.pdf"
  ),
  height = 4,
  width = 4
)
purrr::iwalk(excluded_included_sum, function(df, name) {
  if (base::unique(df$platform) %in% c("droplet", "plate_barcode")) {
    plot <- df |>
      ggplot2::ggplot(aes(x = ic_id_sample, y = n, fill = excluded)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = c("orange", "blue")) +
      ggplot2::labs(
        y = "Number of barcodes",
        fill = "Status",
        title = name
      ) +
      my_theme() +
      ggplot2::theme(
        aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(0.9, 1),
        legend.key.size = unit(0.1, "inch"),
        legend.spacing.y = unit(0, "inch"),
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(margin = margin(l = 0)),
        legend.title = element_blank()
      )
  } else if (base::unique(df$platform) %in% c("plate")) {
    plot <- df |>
      ggplot2::ggplot(aes(x = ic_id_donor, y = n, fill = excluded)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = c("orange", "blue")) +
      ggplot2::labs(
        y = "Number of barcodes",
        fill = "Status",
        title = name
      ) +
      my_theme() +
      ggplot2::theme(
        aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(0.9, 1),
        legend.key.size = unit(0.1, "inch"),
        legend.spacing.y = unit(0, "inch"),
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(margin = margin(l = 0)),
        legend.title = element_blank()
      )
  }
  print(plot)
})
dev.off()