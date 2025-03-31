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
# ## Check if any samples have NA values in quality metrics
# na_values <- quality_met |> modify_depth(
#   1,
#   ~ dplyr::select(., barcode, qc_metrics) |> 
#     dplyr::filter(if_any(c(-barcode), is.na)) |> 
#     pivot_longer(c(-barcode), names_to = "column", values_to = "logical")) |> 
#   bind_rows(.id = "study") |> 
#   pivot_wider(id_cols = c(barcode, study), names_from = column, values_from = logical)
# 
# quality_met[unique(na_values$study)] |> 
#   bind_rows() |> 
#   dplyr::filter(barcode %in% na_values$barcode) |> 
#   dplyr::select(barcode, "Uniquely_mapped_reads_%", "%_of_reads_unmapped_other", "%_of_reads_unmapped_too_short", "%_of_reads_unmapped_too_many_mismatches") |> 
#   tidyr::pivot_longer(cols = -c(barcode), names_to = "qc", values_to = "values") |> 
#   dplyr::mutate(qc = factor(qc, levels = c("Uniquely_mapped_reads_%", "%_of_reads_unmapped_other", "%_of_reads_unmapped_too_short", "%_of_reads_unmapped_too_many_mismatches"))) |> 
#   ggplot2::ggplot(aes(y = barcode, x = values)) +
#   ggplot2::geom_bar(stat = "identity", position = position_dodge()) +
#   facet_wrap(~qc, nrow = 1, scales = "free_x")+
#   my_theme() 

purrr::imap(quality_met, function(df, name) {
  pdf(file = base::paste0(here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/"), name, "_quality_metrics.pdf"), height = 2, width = 2)
  df |> 
    dplyr::select(dplyr::all_of(qc_metrics)) |> 
    purrr::imap(~ print(plot_hist_qc(.x, .y)))
  dev.off()
  })

purrr::imap(quality_met, function(df, name) {
  pdf(file = base::paste0(here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/"), name, "_quality_metrics_per_sample.pdf"), height = 2, width = 2)
  df |> 
    dplyr::select(dplyr::all_of(qc_metrics)) |> 
    purrr::imap(~ print(plot_hist_qc(.x, .y)))
  dev.off()
})


