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

quality_met <- purrr::map(paths, ~vroom::vroom(.x, delim = ",", col_names = TRUE))

# Quality metrics ---------------------------------------------------------

purrr::imap(quality_met, function(df, name) {
  pdf(file = base::paste0(here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/"), name, "_quality_metrics.pdf"), height = 2, width = 2)
  if(base::unique(df$platform) %in% c("droplet", "plate_barcode")){
    df |>
      dplyr::select(dplyr::all_of(qc_metrics_droplet)) |>
      purrr::imap(~ print(plot_hist_qc(.x, .y)))
  } else {
    df |>
      dplyr::select(dplyr::all_of(qc_metrics_plate)) |>
      purrr::imap(~ print(plot_hist_qc(.x, .y)))
  }
  dev.off()
})

## Star quality
purrr::imap(quality_met, function(df, name) {
  pdf(file = base::paste0(here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/"), name, "_quality_star.pdf"), height = 2, width = 2)
  df |> 
    dplyr::rowwise() %>%
    dplyr::mutate("Unmapped_reads_%" = sum(dplyr::c_across(tidyselect::starts_with("%_of_reads_unmapped_")), na.rm = TRUE)) %>%
    dplyr::ungroup()|>  
    dplyr::select(ic_id, dplyr::all_of(qc_star)) |> 
    dplyr::distinct() |> 
    dplyr::select(dplyr::all_of(qc_star)) |> 
    purrr::imap(~ print(plot_hist_star(.x, .y)))
  dev.off()
})

