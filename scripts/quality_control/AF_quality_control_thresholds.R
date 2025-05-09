# Description -------------------------------------------------------------
# From quality control plots, add threshold data for each study and specific donors / samples
# We will also add thresholds and annotate whether a sample or donor is excluded

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
                        # Removed failed samples from barcode
                        dplyr::filter(!ic_id %in%  c("ic_25_11_11", "ic_25_7_7"), !donor == "excluded"))
## Thresholds
# Split into a list
thresholds <- vroom::vroom(here::here("islet_cartography_scrna/data/quality_control/first_pass/first_pass_threshold.csv"), 
                           delim = ";", 
                           col_names = TRUE) |> 
  (function(df) base::split(df, base::factor(df$name)))()

# Not split into a list
thresholds_df <- vroom::vroom(here::here("islet_cartography_scrna/data/quality_control/first_pass/first_pass_threshold.csv"), 
                              delim = ";", 
                              col_names = TRUE)

# Preprocess --------------------------------------------------------------
# Divide Kang into cell and nuclei
quality_met[["Kang_cell"]] <- quality_met[["Kang"]] |>dplyr::filter(cell_nuclei == "cell")
quality_met[["Kang_nuclei"]] <- quality_met[["Kang"]] |>dplyr::filter(cell_nuclei == "nuclei")
quality_met[["Kang"]] <- NULL

# Add column which states is a donor / sample is excluded -----------------
quality_met <- quality_met |> 
  purrr::modify_depth(1, ~ dplyr::mutate(.,
                                         "Excluded" = dplyr::case_when(ic_id %in% c("ic_25_11_11", # poor rank plot
                                                                                    "ic_25_7_7", # poor rank plot
                                                                                    "ic_11_12_359" # Only one cell
                                                                                    ) ~ "excluded",
                                                                       donor == "excluded" ~ "excluded",
                                                                       .default = "included")))

# Add donor specific thresholds -------------------------------------------
sample_thresholds <- data.frame("name" = c(rep("Gurp",2), rep("HPAP_10x", 4), rep("Wang_Sander", 4), rep("Kang_nuclei", 2)),
                                "ic_id_sample" = c("5", "6", "9", "30", "13", "14", "4", "3", "16", "19", "4", "6"),
                                "threshold_nUMIs_lower" = c(10000, 10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1500, 1500),
                                "threshold_nUMIs_upper" = c(900000, 50000, 50000, 15000, 50000, 50000, 10000, 10000, 10000, 10000, 25000, 25000),
                                "threshold_nFeatures_lower" = c(3000, 3000, 1000, 1000, 1000, 2000, 500, 500, 500, 500, 500, 500),
                                "threshold_nFeatures_upper" = c(9000, 8000, 7000, 7000, 7000, 9000, 10000, 10000, 10000, 10000, 10000, 10000),
                                "threshold_complexity_lower" = c(0.3, 0.1, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.4, 0.2, 0.2),
                                "threshold_complexity_upper" = c(0.8, 0.55, 1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4))

donor_thresholds <- data.frame("name" = c("HPAP_patch_22"),
                               "ic_id_sample" = c("5"),
                               "threshold_complexity_lower" = c(0.002),
                               "threshold_complexity_upper" = c(0.015))


