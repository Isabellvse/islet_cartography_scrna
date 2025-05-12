# Description -------------------------------------------------------------
# Create plots for quality control metrics for each study

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
## Quality metrics ----
paths <- base::list.files(path = here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics"),
                          pattern = "quality_metrics.csv",
                          full.names = TRUE, 
                          recursive = TRUE)

base::names(paths) <- purrr::map_chr(paths, ~ stringr::str_extract(.x, pattern = "(?<=quality_metrics/)[^/]+(?=_quality_metrics\\.csv)"))

quality_met <- purrr::map(paths, ~vroom::vroom(.x, delim = ",", col_names = TRUE)) %>% 
  purrr::modify_depth(1, ~dplyr::rowwise(.) %>%
                          dplyr::mutate("Unmapped_reads_%" = sum(dplyr::c_across(tidyselect::starts_with("%_of_reads_unmapped_")), na.rm = TRUE)) %>%
                          dplyr::ungroup() %>% 
                        # Removed failed samples from find empty droplets, so they excluded from this analysis
                        dplyr::filter(!ic_id %in%  c("ic_25_11_11", "ic_25_7_7"), !donor == "excluded"))
## Thresholds
thresholds <- vroom::vroom(here::here("islet_cartography_scrna/data/quality_control/first_pass/first_pass_threshold.csv"), 
                           delim = ";", 
                           col_names = TRUE) |> 
  (function(df) base::split(df, base::factor(df$name)))()

thresholds_df <- vroom::vroom(here::here("islet_cartography_scrna/data/quality_control/first_pass/first_pass_threshold.csv"), 
                           delim = ";", 
                           col_names = TRUE)

# Preprocess --------------------------------------------------------------
# Divide Kang into cell and nuclei
quality_met[["Kang_cell"]] <- quality_met[["Kang"]] |>dplyr::filter(cell_nuclei == "cell")
quality_met[["Kang_nuclei"]] <- quality_met[["Kang"]] |>dplyr::filter(cell_nuclei == "nuclei")
quality_met[["Kang"]] <- NULL


# Plot QC metrics ---------------------------------------------------------
pdf(
  file = here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots.pdf"),
  height = 2,
  width = 2
)
purrr::iwalk(quality_met, function(df, name) {
  subtitle <- name
  
  if (base::unique(df$platform) %in% c("droplet", "plate_barcode")){
    cols <- qc_met_thres_droplet
  } else if (base::unique(df$platform) %in% c("plate")) {
    cols <- qc_met_thres_plate
  }
  
  df |> 
    dplyr::select(tidyselect::all_of(cols)) |>
    purrr::iwalk(~ print(plot_hist_qc_thres(.x, .y, lower = NULL, upper = NULL, subtitle = subtitle)))
  
})
dev.off()

# Zoom on some QC plots ---------------------------------------------------
zoom_list <- base::list("Dai" = base::data.frame("nCounts" = c(0,300000),
                                                 "complexity" = c(0, 0.05)),
                        "Enge" = base::data.frame("nCounts" = c(0, 1000000),
                                                  "complexity" = c(0, 0.02)),
                        "Fang" = base::data.frame("nUMIs" = c(0, 1000),
                                                  "nFeatures" = c(0, 1000)),
                        "Gurp" = base::data.frame("nUMIs" = c(0, 40000)),
                        "HPAP_10x" = base::data.frame("nUMIs" = c(0, 50000),
                                                      "contrast_fraction" = c(0, 1.5)),
                        "HPAP_patch_22" = base::data.frame("nCounts" = c(0, 500000)),
                        "HPAP_patch_23" = base::data.frame("nCounts" = c(0, 600000),
                                                           "complexity" = c(0, 0.01)),
                        "Kang_cell" = base::data.frame("nUMIs" = c(0, 30000)),
                        "Kang_nuclei" = base::data.frame("nUMIs" = c(0, 10000)),
                        "Lawlor" = base::data.frame("complexity" = c(0, 0.02)),
                        "Mauvais_Jarvis" = base::data.frame("nUMIs" = c(0, 40000),
                                                            "contrast_fraction" = c(0, 2)),
                        "Motakis" = base::data.frame("nUMIs" = c(0, 10000)),
                        "Segerstolpe" = base::data.frame("nCounts" = c(0, 500000)),
                        "Shrestha" = base::data.frame("nUMIs" = c(0, 10000)),
                        "Wang_Sander" = base::data.frame("nUMIs" = c(0, 10000),
                                                         "nFeatures" = c(0, 5000),
                                                         "contrast_fraction" =c(0, 0.5)),
                        "Xin_Diabetes" = base::data.frame("nUMIs" = c(0, 10000)),
                        "Xin" = base::data.frame("nCounts" = c(0, 750000)),
                        "Zhang" = base::data.frame("nUMIs" = c(0, 1500),
                                                   "nFeatures" = c(0, 500)),
                        "Camunas" = base::data.frame("nCounts" = c(0, 600000),
                                                     "nFeatures" = c(0, 4000)))

pdf(
  file = here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_zoomed.pdf"),
  height = 2,
  width = 2
)

purrr::iwalk(quality_met, function(df, name) {
  subtitle <- name
  zooms <- zoom_list[[name]]
  
  # skip if there's no zoom info
  if (base::is.null(zooms)) return()
  
  # plot only the variables in zoom_list for this dataset
  intersecting_vars <- base::intersect(base::names(zooms), base::colnames(df))
  
  purrr::walk(intersecting_vars, function(varname) {
    xlim <- base::as.numeric(zooms[[varname]])
    plot_data <- df[[varname]]
    
    print(plot_hist_qc_thres(
      .x = plot_data,
      .y = varname,
      lower = NULL,
      upper = NULL,
      subtitle = subtitle,
      .xlim = xlim
    ))
  })
})

dev.off()

# Plot with thresholds ----------------------------------------------------
## Passed or not passed ----
# 0 = passed 
# 1 = failed
checked <- purrr::imap(quality_met, function(qc, name) {
  check_thresholds(qc_met = qc, thresholds = thresholds[[name]])
})


## Plot ----
pdf(
  file = here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_thresholds.pdf"),
  height = 2,
  width = 2
)
purrr::iwalk(quality_met, function(df, name) {
  
  thresholds_filtered <- thresholds[[name]]
  subtitle <- name
  
  # Number of barcodes that fail quality control
  df_thres <- checked[[name]] |> 
    dplyr::select(tidyselect::any_of(base::paste0(c(qc_met_thres_plate, qc_met_thres_droplet), "_thres"))) |> 
    dplyr::summarise_all(sum) |> 
    dplyr::rename_all(~stringr::str_replace(.,"_thres","")) 
  
  # Number of barcodes that fail more than one quality control
  df_multi <- checked[[name]] |> 
    dplyr::select(dplyr::ends_with("thres_multipass")) |> 
    dplyr::summarise_all(sum) |> 
    dplyr::rename_all(~stringr::str_replace(.,"_thres_multipass","")) 
  
  # Total number of barcodes
  n_total <- checked[[name]] |> base::nrow()
  
  if (base::unique(df$platform) %in% c("droplet", "plate_barcode")){
    cols <- qc_met_thres_droplet
  } else if (base::unique(df$platform) %in% c("plate")) {
    cols <- qc_met_thres_plate
  }

    df |>
      dplyr::select(tidyselect::all_of(cols)) |>
      purrr::iwalk(~ {
        lower_col <- paste0("threshold_", .y, "_lower")
        upper_col <- paste0("threshold_", .y, "_upper")
        lower <- thresholds_filtered[[lower_col]][1]
        upper <- thresholds_filtered[[upper_col]][1]
        
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
        caption <- base::paste0("Total: ", n_total, ", Failed: ", fail, ", Failed multi: ", round((fail_multi/fail)*100,2), "%")
        
        # Plot
        base::print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle, caption = caption))
        }
        )
    }
  )
dev.off()


# Barcodes removed per sample / donor -------------------------------------
# Plot number of barcodes removed from each sample (droplet based) or donor (plate based)

## Plot ----
pdf(
  file = here::here(
    "islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_per_donor_sample_bar.pdf"
  ),
  height = 3,
  width = 16  # Adjust width to accommodate multiple facets in one row
)
purrr::iwalk(checked, function(df, subtitle) {
  subtitle = subtitle
  
  if (base::unique(df$platform) %in% c("droplet", "plate_barcode")) {
    df_thres <- df |>
      dplyr::select(ic_id_sample, tidyselect::all_of(base::paste0(qc_met_thres_droplet, "_thres"))) |>
      dplyr::group_by(ic_id_sample) |>
      dplyr::summarise_all(sum) |>
      dplyr::rename_all( ~ stringr::str_replace(., "_thres", ""))
    
    # Total number of barcodes
    n_total <- df |>
      dplyr::group_by(ic_id_sample) |>
      dplyr::tally() |>
      dplyr::full_join(df_thres, by = "ic_id_sample") |>
      dplyr::mutate(dplyr::across(
        c(-ic_id_sample, -n),
        .fns = function(qc) {
          round((qc / n) * 100, 2)
        }
      ),
      ic_id_sample = as.character(ic_id_sample)) |>
      dplyr::select(-n) |>
      tidyr::pivot_longer(-ic_id_sample, names_to = "qc", values_to = "value")
    
    # Identify rows with missing values
    missing_values <- n_total  |>  dplyr::filter(is.na(value))
    print("Study:")
    print(subtitle)
    print("Rows with missing values:")
    print(missing_values)
    
    # Identify rows with values outside the scale range (assuming scale range is 0-100)
    outside_scale_range <- n_total  |>  dplyr::filter(value < 0 | value > 100)
    print("Rows with values outside the scale range:")
    print(outside_scale_range)
    
    plot <- n_total |>
      ggplot2::ggplot(aes(y = ic_id_sample, x = value)) +
      ggplot2::geom_bar(stat = "identity",
                        position = position_dodge(),
                        fill = "black") +
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
    
    n_total_wide <- n_total %>% tidyr::pivot_wider(id_cols = ic_id_sample, names_from = qc, values_from = value)
    print(n_total_wide)
  } else if (base::unique(df$platform) %in% c("plate")) {
    df_thres <- df |>
      dplyr::select(ic_id_donor, tidyselect::all_of(base::paste0(qc_met_thres_plate, "_thres"))) |>
      dplyr::group_by(ic_id_donor) |>
      dplyr::summarise_all(sum) |>
      dplyr::rename_all( ~ stringr::str_replace(., "_thres", ""))
    
    # Total number of barcodes
    n_total <- df |>
      dplyr::group_by(ic_id_donor) |>
      dplyr::tally() |>
      dplyr::full_join(df_thres, by = "ic_id_donor") |>
      dplyr::mutate(dplyr::across(
        c(-ic_id_donor, -n),
        .fns = function(qc) {
          round((qc / n) * 100, 2)
        }
      ),
      ic_id_donor = as.character(ic_id_donor)) |>
      dplyr::select(-n) |>
      tidyr::pivot_longer(-ic_id_donor,
                          names_to = "qc",
                          values_to = "value")
    
    # Identify rows with missing values
    missing_values <- n_total |>  dplyr::filter(is.na(value))
    print("Study:")
    print(subtitle)
    print("Rows with missing values:")
    print(missing_values)
    
    # Identify rows with values outside the scale range (assuming scale range is 0-100)
    outside_scale_range <- n_total  |>  dplyr::filter(value < 0 | value > 100)
    print("Rows with values outside the scale range:")
    print(outside_scale_range)
    
    plot <- n_total |>
      ggplot2::ggplot(aes(y = ic_id_donor, x = value)) +
      ggplot2::geom_bar(stat = "identity",
                        position = position_dodge(),
                        fill = "black") +
      ggplot2::geom_text(aes(label = value), size = 1.5, hjust = -0.2) +
      ggplot2::labs(
        subtitle = subtitle,
        x = "% of total barcodes",
        y = "Donor"
      ) +
      ggplot2::scale_x_continuous(limits = c(0, 100), breaks = scales::pretty_breaks(n = 5)) +
      ggplot2::facet_wrap( dplyr::case_when(
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
        TRUE ~ qc) ~ ., nrow = 1) + 
      my_theme()
    
    print(plot)
    
    n_total_wide <- n_total %>% tidyr::pivot_wider(id_cols = ic_id_donor, names_from = qc, values_from = value)
    print(n_total_wide)
  }
}
)
dev.off()


# QC plots per donor / sample with thresholds -----------------------------

# Split qc values by donor or sample according to platform
qc_split <- quality_metqc_split <- quality_met |> 
  purrr::modify_depth(1, ~ if (base::unique(.x$platform) %in% c("droplet", "plate_barcode")) {
    base::split(.x, .x$ic_id_sample)
  } else if (base::unique(.x$platform) %in% c("plate")) {
    base::split(.x, .x$ic_id_donor)
}) |> 
  purrr::list_flatten() %>% 
  purrr::modify_at(c("Kang_cell_1", "Kang_cell_3", "Kang_cell_5"), ~dplyr::mutate(., name = paste0(name, "_cell"))) %>% 
  purrr::modify_at(c("Kang_nuclei_2", "Kang_nuclei_4", "Kang_nuclei_6"), ~dplyr::mutate(., name = paste0(name, "_nuclei")))


# Add thresholds values to qc values by donor / sample
thresholds_split <- qc_split |> 
  purrr::modify_depth(1, ~ 
                        if (base::unique(.x$platform) %in% c("droplet", "plate_barcode")) {
                          dplyr::select(.x, id = ic_id_sample, name) |> 
                            dplyr::distinct()
                        } else if (base::unique(.x$platform) %in% c("plate")) {
                          dplyr::select(.x, id = ic_id_donor, name) |> 
                            dplyr::distinct()
                        }) |> 
  dplyr::bind_rows(.id = "id") |> 
  dplyr::left_join(y = thresholds_df, by = "name")

## Plot ----
pdf(
  file = here::here(
    "islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_per_donor_sample.pdf"),
  height = 2,
  width = 2
)
purrr::iwalk(qc_split, function(df, name) {
  thresholds_filtered <- dplyr::filter(thresholds_split, id == !!name)
  subtitle <- name
  if (base::unique(df$platform) %in% c("droplet", "plate_barcode")) {
    df |>
      dplyr::select(dplyr::all_of(qc_met_thres_droplet)) |>
      purrr::imap(~ {
        lower_col <- paste0("threshold_", .y, "_lower")
        upper_col <- paste0("threshold_", .y, "_upper")
        lower <- thresholds_filtered[[lower_col]][1]
        upper <- thresholds_filtered[[upper_col]][1]
        
        # Check if lower or upper is NA or NULL and replace accordingly
        if (is.na(lower) || is.null(lower)) {
          lower <- NULL
        }
        
        if (is.na(upper) || is.null(upper)) {
          upper <- NULL
        }
        
        print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle))
      })
  } else if (base::unique(df$platform) %in% c("plate")) {
    df |>
      dplyr::select(dplyr::all_of(qc_met_thres_plate)) |>
      purrr::imap(~ {
        lower_col <- paste0("threshold_", .y, "_lower")
        upper_col <- paste0("threshold_", .y, "_upper")
        lower <- thresholds_filtered[[lower_col]][1]
        upper <- thresholds_filtered[[upper_col]][1]
        
        # Check if lower or upper is NA or NULL and replace accordingly
        if (is.na(lower) || is.null(lower)) {
          lower <- NULL
        }
        
        if (is.na(upper) || is.null(upper)) {
          upper <- NULL
        }
        
        print(plot_hist_qc_thres(.x, .y, lower, upper, subtitle = subtitle))
      })
  }
})
dev.off()




