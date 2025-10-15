# Description -------------------------------------------------------------
# Here we do quality control per leiden cluster, in order to remove low quality clusters and doublets

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
adata_obs <- vroom::vroom(here::here("islet_cartography_scrna/data/integrate/second_pass/files/adata_obs.csv"))

# Preprocess --------------------------------------------------------------
# Only keep columns of interest
adata_flt <- adata_obs |>  dplyr::select(n_umi, n_feature, mitochondrial_fraction, n_count, 
                          tidyselect::starts_with("azimuth"), doublet_probability, leiden_res_5.50)

# Calculate mean and median values per leiden cluster
adata_sum <- adata_flt |> 
  dplyr::group_by(leiden_res_5.50) |> 
  dplyr::summarise(dplyr::across(base::is.numeric,
                                 list(mean = ~mean(.x, na.rm = TRUE),
                                      median = ~median(.x, na.rm = TRUE))),
                   n_cells = n(),
                   .groups = "drop")


# Distances in gene set scores --------------------------------------------
# Identify which leiden cluster that has the max score
max_score_col <- adata_sum %>%
  dplyr::select(leiden_res_5.50, tidyselect::matches("^azimuth.*_score_median$")) %>%
  tidyr::pivot_longer(cols = -leiden_res_5.50, names_to = "score_col", values_to = "value") %>%
  dplyr::group_by(leiden_res_5.50) %>%
  dplyr::slice_max(order_by = value, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(leiden_res_5.50, score_col)

# Join max score info back to adata_sum
adata_sum_max <- adata_sum %>%
  dplyr::left_join(max_score_col, by = "leiden_res_5.50")

# Calculate absolute differences from the max score
adata_sum_dist <- adata_sum_max %>%
  dplyr::rowwise() %>%
  dplyr::mutate(across(
    tidyselect::matches("^azimuth.*_score_median$"),
    ~ abs(. - get(score_col)),
    .names = "{.col}_dist_from_max"
  )) %>%
  dplyr::ungroup()

# Combine per barcode, and calculate perc of cells per cluster that have a doublet probability higher than 70 %
test <- adata_sum_dist |>
  dplyr::select(leiden_res_5.50, tidyselect::contains("_dist_from_max"), score_col) |> 
  tidyr::pivot_longer(c(-leiden_res_5.50, -score_col, -score_col), ) |> 
  dplyr::group_by(leiden_res_5.50) |> 
  dplyr::mutate(remove = case_when(value == 0 & score_col == stringr::str_remove(name, "_dist_from_max") ~ TRUE,
                                   stringr::str_remove(name, "_dist_from_max") == "azimuth_cycling_score_median" ~ TRUE,
                                   .default = FALSE)) |> 
  dplyr::filter(!remove == TRUE) |> 
  dplyr::select(-remove) |> 
  dplyr::slice_min(order_by = value, n= 1) |> 
  dplyr::ungroup() |> 
  dplyr::right_join(adata_obs |> dplyr::select(leiden_res_5.50, doublet_probability)) |> 
  dplyr::mutate(doublet = dplyr::case_when(doublet_probability >= 0.7 ~ "doublet",
                                           doublet_probability < 0.7 ~ "singlet")) |> 
  dplyr::group_by(leiden_res_5.50) |> 
  dplyr::mutate(perc_doublet = sum(doublet == "doublet", na.rm = TRUE) / dplyr::n() * 100) |> 
  dplyr::ungroup()

# Extract cluster based metrics
test_2 <- test |> 
  dplyr::select(-doublet_probability, -doublet) |> dplyr::distinct()

# Plot perc of doublets vs minimum distance in gene module scores
plot(test_2$perc_doublet, test_2$value)
abline(h=0.8, v =25, col = "red")

# Try to set a threshold, 
# those clusters with a distance less than 0.8 and a doublet % higher than 20
# are marked as doublet clusters
doublet_clusters <- test_2 |> 
  dplyr::filter(value <= 0.8 & perc_doublet >= 25) |> 
  pull(leiden_res_5.50)


adata_flt |> 
  dplyr::mutate(doublet = case_when(leiden_res_5.50 %in% doublet_clusters ~ "doublet",
                                    .default = "singlet")) |> 
  ggplot2::ggplot(aes(x = doublet, y = mitochondrial_fraction)) +
  ggplot2::geom_violin() +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.2, fill = "transparent") 

# QC plots ----------------------------------------------------------------
plot_hist_qc <- function(.x, .y, title) {
  ggplot2::ggplot(data = tibble::tibble(value = .x), ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 100, fill = "black") +
    ggplot2::labs(title = title,
      x = title, 
      y = "Frequency") +
    ggplot2::scale_x_continuous(
      labels = function(x) ifelse(x == 0, x, ifelse(abs(x) >= 10000 | abs(x) < 0.0001, scales::scientific(x, digits = 1), x)),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    my_theme() +
    ggplot2::theme(aspect.ratio = 1)
}

pdf(
  file = here::here("islet_cartography_scrna/data/integrate/second_pass/plot/qc_plots_per_cluster.pdf"),
  height = 2,
  width = 2
)
adata_sum_full |> 
  dplyr::select(tidyselect::contains("_vs_"), tidyselect::contains("_median"), -tidyselect::contains("azimuth")) |>
  purrr::iwalk(~ print(plot_hist_qc(.x, .y, title = .y)))
dev.off()

