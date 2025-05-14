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
quality_met[["Kang_cell"]] <- quality_met[["Kang"]] |>dplyr::filter(cell_nuclei == "cell") |> 
  dplyr::mutate(name = paste0(name, "_cell"))
quality_met[["Kang_nuclei"]] <- quality_met[["Kang"]] |>dplyr::filter(cell_nuclei == "nuclei") |> 
  dplyr::mutate(name = paste0(name, "_nuclei"))
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


# Plots -------------------------------------------------------------------
## Passed or not passed ----
# 0 = passed 
# 1 = failed
checked <- purrr::imap(quality_met, function(qc, name) {
  check_thresholds(qc_met = qc, thresholds = thresholds[[name]])
})

## Histogram per study with study level thresholds ----
pdf(
  file = here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_thresholds.pdf"),
  height = 2,
  width = 2
)
iwalk(quality_met, function(qc_met, name){
  visualize_qc_hist_with_thresholds(qc_met = qc_met, thresholds = thresholds[[name]], checked = checked[[name]], subtitle = name)
  return(checked)
})
dev.off()

## Barcodes removed per sample / donor -------------------------------------
# Plot number of barcodes removed from each sample (droplet based) or donor (plate based)
pdf(
  file = here::here(
    "islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_per_donor_sample_bar.pdf"
  ),
  height = 3,
  width = 16  # Adjust width to accommodate multiple facets in one row
)
purrr::iwalk(checked, function(df, subtitle) {
  subtitle = subtitle
  
  visualize_barcodes_failed_qc(checked = df, subtitle = subtitle)
})
dev.off()

# Combine thresholds, passing and quality control together ----------------
quality_met <- imap(quality_met, function(df, name){
  df |> 
    dplyr::left_join(y = thresholds[[name]], 
                     relationship = "many-to-one") |> 
    dplyr::full_join(y = checked[[name]])
})

## qc metrics per sample / donor ----
pdf(
  file = here::here(
    "islet_cartography_scrna/data/quality_control/first_pass/plots/qc_plots_per_donor_sample.pdf"),
  height = 2,
  width = 2
)

# Split by sample or donor
quality_met_split <- quality_met |> 
  purrr::modify_depth(1, ~ if (base::unique(.x$platform) %in% c("droplet", "plate_barcode")) {
    base::split(.x, .x$ic_id_sample)
  } else if (base::unique(.x$platform) %in% c("plate")) {
    base::split(.x, .x$ic_id_donor)}) |> 
  purrr::list_flatten()

# make hitogram plots for each sample / donor
purrr::iwalk(quality_met_split, function(df, name){
  subtitle <- name
  print(subtitle)
  visualize_qc_hist_threshold_per_sample(qc_met_split = df, subtitle = subtitle)
})

dev.off()
