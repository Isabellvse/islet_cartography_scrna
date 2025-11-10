# Description -------------------------------------------------------------
# Here we do quality control per leiden cluster, in order to remove low quality clusters and doublets

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
adata_obs <- vroom::vroom(here::here("islet_cartography_scrna/data/integrate/second_pass/files/adata_obs.csv"))

# Helper function ----------------------------------------------------------
plot_hist_qc <- function(.x, .y, title) {
  ggplot2::ggplot(data = tibble::tibble(value = .x), ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 50, fill = "black") +
    ggplot2::labs(title = title,
                  x = title, 
                  y = "Frequency") +
    ggplot2::scale_x_continuous(
      labels = function(x) ifelse(x == 0, x, ifelse(abs(x) >= 10000 | abs(x) < 0.0001, scales::scientific(x, digits = 2), x)),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    my_theme() +
    ggplot2::theme(aspect.ratio = 1)
}

# Per cluster analysis ----------------------------------------------------
## Preprocess ----
# Only keep columns of interest
adata_flt <- adata_obs |>  dplyr::select(n_umi, n_feature, mitochondrial_fraction, n_count, coding_fraction, contrast_fraction, complexity,
                                         tidyselect::starts_with("azimuth"), doublet_probability, leiden_res_10.00)

# Calculate mean and median values per leiden cluster
adata_sum <- adata_flt |> 
  dplyr::group_by(leiden_res_10.00) |> 
  dplyr::summarise(dplyr::across(base::is.numeric,
                                 list(mean = ~mean(.x, na.rm = TRUE),
                                      median = ~median(.x, na.rm = TRUE))),
                   n_cells = n(),
                   .groups = "drop")

# Calculate distance in module score --------------------------------------
# Identify which leiden cluster that has the max score
max_score_col <- adata_sum %>%
  dplyr::select(leiden_res_10.00, tidyselect::matches("^azimuth.*_score_median$")) %>%
  tidyr::pivot_longer(cols = -leiden_res_10.00, names_to = "max_score_col", values_to = "value") %>%
  dplyr::group_by(leiden_res_10.00) %>%
  dplyr::slice_max(order_by = value, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(leiden_res_10.00, max_score_col)

# Join max score info back to adata_sum
adata_sum_max <- adata_sum %>%
  dplyr::left_join(max_score_col, by = "leiden_res_10.00")

# Calculate absolute differences from the max score
adata_sum_dist <- adata_sum_max %>%
  dplyr::rowwise() %>%
  dplyr::mutate(across(
    tidyselect::matches("^azimuth.*_score_median$"),
    ~ base::abs(. - base::get(max_score_col)),
    .names = "{.col}_dist_from_max"
  )) %>%
  dplyr::ungroup()

# Get scores between highest and second higest module score and combin with adata_sum
adata_comb <- adata_sum_dist |>
  dplyr::select(leiden_res_10.00, tidyselect::contains("_dist_from_max"), max_score_col) |> 
  tidyr::pivot_longer(c(-leiden_res_10.00, -max_score_col), names_to = "sec_score_col", values_to = "module_difference") |> 
  dplyr::group_by(leiden_res_10.00) |> 
  dplyr::mutate(remove = case_when(module_difference == 0 & max_score_col == stringr::str_remove(sec_score_col, "_dist_from_max") ~ TRUE,
                                   stringr::str_remove(sec_score_col, "_dist_from_max") == "azimuth_cycling_score_median" ~ TRUE,
                                   .default = FALSE)) |> 
  dplyr::filter(!remove == TRUE) |> 
  dplyr::select(-remove) |> 
  dplyr::slice_min(order_by = module_difference, n= 1) |> 
  dplyr::ungroup() |> 
  dplyr::right_join(y = adata_sum)

## Plot module score vs db score -------------------------------------------

# Plot clusters arranged by module score, and sized by their doublet score, colors by if 
# they are defined as a doublet cluster
module_doublet <- adata_comb |> 
  dplyr::mutate(filter = dplyr::case_when(doublet_probability_median >= 0.4 ~ "doublet",
                                          .default = "singlet"),
                only_04 = dplyr::case_when(leiden_res_10.00 %in% diff_clusters ~ "only in 0.4",
                                           .default = "in 0.4 and 0.5")) |> 
  ggplot2::ggplot(aes(x = doublet_probability_median, y = module_difference)) +
  ggplot2::geom_point(aes(color = only_04, shape = filter)) +
  ggplot2::scale_shape_manual(values=c(4, 16)) +
  ggplot2::labs(title = "Module score vs\ndoublet probability",
                subtitle = "Leiden resolution 10.00",
                y = "Median module score difference",
                x = "Median doublet probability",
                color = "Median doublet\nprobability",
                shape = "Doublet classification\nscore >= 0.4") +
  ggplot2::geom_text(aes(label=base::ifelse(only_04 == "only in 0.4", 
                                            as.character(leiden_res_10.00),'')),
                     hjust=0,vjust=0, size = 1) +
  my_theme() +
  theme(aspect.ratio = 1)

pdf(
  file = here::here("islet_cartography_scrna/data/integrate/second_pass/plot/module_scor_vs_doublet_per_cluster_leiden_10.pdf"),
  height = 2,
  width = 2.5
)
module_doublet
dev.off()



# Apply filter ------------------------------------------------------------
# Filter per cluster
filter_cluster <- adata_sum  |>  
  dplyr::mutate(filter = dplyr::case_when(doublet_probability_median >= 0.4 ~ 1,
                                          n_umi_median <= 2000 ~ 1,
                                          contrast_fraction_median <= 0.4 ~ 1,
                                          complexity_median <= 0.1 ~ 1,
                                          .default = 0)) |> 
  dplyr::filter(filter >= 1) |> 
  dplyr::pull(leiden_res_10.00)



# Filter per cell
filter_cells <- adata_obs |> 
  dplyr::mutate(filter = dplyr::case_when(doublet_probability >= 0.9 ~ 1,
                                          leiden_res_10.00 %in% filter_cluster ~ 1,
                                          .default = 0)) |> 
  dplyr::filter(filter < 1) |> 
  dplyr::select(barcode)

cells_remove <- adata_obs |> 
  dplyr::mutate(filter = dplyr::case_when(doublet_probability >= 0.9 ~ 1,
                                          leiden_res_10.00 %in% filter_cluster ~ 1,
                                          .default = 0)) |> 
  dplyr::filter(filter >= 1) |> 
  dplyr::select(barcode)

# sum of cells in each cluster following qc
sum_of_cells <- adata_obs |> 
  dplyr::select(barcode, leiden_res_10.00) |> 
  dplyr::filter(barcode %in% filter_cells$barcode) |> 
  dplyr::group_by(leiden_res_10.00) |> 
  dplyr::summarise(
    n_cells = n(),
    .groups = "drop")

sum_of_cells_anno <- adata_obs |> 
  dplyr::select(barcode, study_cell_annotation_harmonized) |> 
  dplyr::filter(barcode %in% filter_cells$barcode) |> 
  dplyr::group_by(study_cell_annotation_harmonized) |> 
  dplyr::summarise(
    n_cells = n(),
    .groups = "drop")

# Save barcodes to be removed ---------------------------------------------
vroom::vroom_write(filter_cells, 
                   here::here("islet_cartography_scrna/data/integrate/second_pass/files/barcodes_keep.csv"))

vroom::vroom_write(cells_remove, 
                   here::here("islet_cartography_scrna/data/integrate/second_pass/files/barcodes_keep.csv"))

## QC plots ----------------------------------------------------------------
pdf(
  file = here::here("islet_cartography_scrna/data/integrate/second_pass/plot/qc_plots_per_cluster_leiden_10.pdf"),
  height = 2,
  width = 2
)
adata_sum |> 
  dplyr::select(tidyselect::contains("_median"), -tidyselect::contains("azimuth")) |>
  purrr::iwalk(~ print(plot_hist_qc(.x, .y, title = .y)))
dev.off()



# QC plots with thresholds ------------------------------------------------
plot_list <-adata_sum %>% 
  dplyr::select(tidyselect::contains("_median"), -tidyselect::contains("azimuth")) |>
  purrr::imap(~ plot_hist_qc(.x, .y, title = .y))

plot_list_2 <- plot_list  |>  
  purrr::imap(function(.x, .y){
    vec <- adata_sum  |>  dplyr::filter(contrast_fraction_median < 0.7)  |>  dplyr::pull(.y) 
    output <- .x + ggplot2::geom_vline(xintercept = vec, color = "red") 
    return(output)
  })

wrap_plots(plot_list_2, nrow = 1) + plot_annotation(title = "Resolution = 10.00")

pdf(
  file = here::here("islet_cartography_scrna/data/integrate/second_pass/plot/qc_plots_per_cluster_leiden_10_with_threshold.pdf"),
  height = 2,
  width = 2
)
# This is lower thresholds
plot_list[[1]] + ggplot2::geom_vline(xintercept = 2000, color = "red") 
plot_list[[6]] + ggplot2::geom_vline(xintercept = 0.4, color = "red") 
plot_list[[7]] + ggplot2::geom_vline(xintercept = 0.05, color = "red") 
plot_list[[8]] + ggplot2::geom_vline(xintercept = 0.4, color = "red")
dev.off()

# Histogram of cells removed with cluster doublet prop filer --------------
doublet_clusters <- adata_sum |> 
  dplyr::filter(doublet_probability_median >= 0.4) |> 
  dplyr::pull(leiden_res_10.00)

pdf(
  file = here::here("islet_cartography_scrna/data/integrate/second_pass/plot/qc_plots_per_cell_doublet_probability.pdf"),
  height = 2,
  width = 2
)
adata_obs |> 
  dplyr::filter(leiden_res_10.00 %in% doublet_clusters) |> 
  ggplot2::ggplot(aes(x = doublet_probability)) + 
  ggplot2::geom_histogram(bins = 50, fill = "black") +
  ggplot2::geom_boxplot(aes(y = -100, group = 1),
               width = 100, 
               outlier.size = 0.5) +
  ggplot2::labs(x = "Doublet probability",
                y = "Frequency",
                title = "Per cell doublet probability") +
  my_theme()
dev.off()

pdf(
  file = here::here("islet_cartography_scrna/data/integrate/second_pass/plot/qc_plots_per_cell_doublet_probability_threshold.pdf"),
  height = 2,
  width = 2
)
# This is lower thresholds
adata_obs |> 
  dplyr::filter(leiden_res_10.00 %in% doublet_clusters) |> 
  ggplot2::ggplot(aes(x = doublet_probability)) + 
  ggplot2::geom_histogram(bins = 50, fill = "black") +
  ggplot2::geom_boxplot(aes(y = -100, group = 1),
                        width = 100, 
                        outlier.size = 0.5) +
  ggplot2::geom_vline(xintercept = 0.9, color = "blue") +
  ggplot2::labs(x = "Doublet probability",
                y = "Frequency",
                title = "Per cell doublet probability") +
  my_theme()
dev.off()

adata_obs |> 
  dplyr::filter(leiden_res_10.00 %in% doublet_clusters) |> 
  dplyr::pull(doublet_probability) |> 
  stats::quantile(na.rm = TRUE)



# Define cells to remove --------------------------------------------------
# Filter per cluster
filter_cluster <- adata_sum  |>  
  dplyr::mutate(filter = dplyr::case_when(doublet_probability_median >= 0.4 ~ 1,
                                          n_umi_median <= 2000 ~ 1,
                                          contrast_fraction_median <= 0.4 ~ 1,
                                          complexity_median <= 0.1 ~ 1,
                                          .default = 0)) |> 
  dplyr::filter(filter >= 1) |> 
  dplyr::pull(leiden_res_10.00)



# Filter per cell
filter_cells <- adata_obs |> 
  dplyr::mutate(filter = dplyr::case_when(doublet_probability >= 0.9 ~ 1,
                                          leiden_res_10.00 %in% filter_cluster ~ 1,
                                          .default = 0)) |> 
  dplyr::filter(filter < 1) |> 
  dplyr::select(barcode)

cells_remove <- adata_obs |> 
  dplyr::mutate(filter = dplyr::case_when(doublet_probability >= 0.9 ~ 1,
                                          leiden_res_10.00 %in% filter_cluster ~ 1,
                                          .default = 0)) |> 
  dplyr::filter(filter >= 1) |> 
  dplyr::select(barcode)

cells_keep <- adata_obs |> 
  dplyr::mutate(filter = dplyr::case_when(doublet_probability >= 0.9 ~ 1,
                                          leiden_res_10.00 %in% filter_cluster ~ 1,
                                          .default = 0)) |> 
  dplyr::filter(filter == 0) |> 
  dplyr::select(barcode)



# Save --------------------------------------------------------------------
# Save csv files
vroom::vroom_write(cells_remove, here::here("islet_cartography_scrna/data/integrate/second_pass/files/barcodes_failed_qc.csv"),
delim = ",",
col_names = TRUE)

vroom::vroom_write(cells_keep, here::here("islet_cartography_scrna/data/integrate/second_pass/files/barcodes_pass_qc.csv"),
                   delim = ",",
                   col_names = TRUE)



# A closer look on clusters in the middle ---------------------------------

pdf(width = 6, height = 2)
adata_comb |>  
  dplyr::filter(leiden_res_10.00 %in% c("53",  "75", "176")) |> 
  dplyr::select(leiden_res_10.00, tidyselect::matches("^azimuth.*_score_median$")) |> 
  dplyr::mutate(leiden_res_10.00 = as.character(leiden_res_10.00))  |>  
  tidyr::pivot_longer(-leiden_res_10.00)  |>  
  dplyr::mutate(split_by = leiden_res_10.00) %>%
  collapse::rsplit(~ split_by) |> 
  purrr::imap(~{
    plot <- .x |> 
      dplyr::mutate(name = stringr::str_remove(name, "azimuth_") |> stringr::str_remove("_score_median"),
                    name = forcats::fct_reorder(name, value, .desc = FALSE)) |> 
      ggplot2::ggplot(aes(y = name, x = value)) +
      ggplot2::geom_bar(stat = "identity", position = position_dodge()) +
      ggplot2::labs(title = .y, 
                    x = "Median module score",
                    y = "Geneset") +
      my_theme() +
      theme(legend.position = "none",
            aspect.ratio = 1)
    return(plot)
  }) |> 
  patchwork::wrap_plots()
dev.off()

## What does the cluster with the highest doublet score look like?
pdf(width = 6, height = 2)
adata_comb |>  
  dplyr::top_n(n = 3, wt = doublet_probability_median) |> 
  dplyr::select(leiden_res_10.00, tidyselect::matches("^azimuth.*_score_median$")) |> 
  dplyr::mutate(leiden_res_10.00 = as.character(leiden_res_10.00))  |>  
  tidyr::pivot_longer(-leiden_res_10.00)  |>  
  dplyr::mutate(split_by = leiden_res_10.00) %>%
  collapse::rsplit(~ split_by) |> 
  purrr::imap(~{
    plot <- .x |> 
      dplyr::mutate(name = stringr::str_remove(name, "azimuth_") |> stringr::str_remove("_score_median"),
                    name = forcats::fct_reorder(name, value, .desc = FALSE)) |> 
      ggplot2::ggplot(aes(y = name, x = value)) +
      ggplot2::geom_bar(stat = "identity", position = position_dodge()) +
      ggplot2::labs(title = .y, 
                    x = "Median module score",
                    y = "Geneset") +
      my_theme() +
      theme(legend.position = "none",
            aspect.ratio = 1)
    return(plot)
  }) |> 
  patchwork::wrap_plots()
dev.off()

## What does the cluster with the lowest doublet score look like?
pdf(width = 6, height = 2)
adata_comb |>  
  dplyr::top_n(n = -3, wt = doublet_probability_median) |> 
  dplyr::select(leiden_res_10.00, tidyselect::matches("^azimuth.*_score_median$")) |> 
  dplyr::mutate(leiden_res_10.00 = as.character(leiden_res_10.00))  |>  
  tidyr::pivot_longer(-leiden_res_10.00)  |>  
  dplyr::mutate(split_by = leiden_res_10.00) %>%
  collapse::rsplit(~ split_by) |> 
  purrr::imap(~{
    plot <- .x |> 
      dplyr::mutate(name = stringr::str_remove(name, "azimuth_") |> stringr::str_remove("_score_median"),
                    name = forcats::fct_reorder(name, value, .desc = FALSE)) |> 
      ggplot2::ggplot(aes(y = name, x = value)) +
      ggplot2::geom_bar(stat = "identity", position = position_dodge()) +
      ggplot2::labs(title = .y, 
                    x = "Median module score",
                    y = "Geneset") +
      my_theme() +
      theme(legend.position = "none",
            aspect.ratio = 1)
    return(plot)
  }) |> 
  patchwork::wrap_plots()
dev.off()

pdf(width = 2, height = 2)
adata_comb |>  
  dplyr::top_n(n = -3, wt = doublet_probability_median)

adata_comb |> dplyr::filter(leiden_res_10.00 %in% c("53",  "75", "176", "137", "138", "156", "92", "136", "162")) |> 
  dplyr::select(leiden_res_10.00, n_umi_median, n_feature_median, n_count_median, doublet_probability_median) |> 
  dplyr::mutate(leiden_res_10.00 = as.character(leiden_res_10.00),
                leiden_res_10.00 = forcats::fct_reorder(leiden_res_10.00, doublet_probability_median, .desc = FALSE)) |> 
  ggplot(aes(y = n_feature_median, x = leiden_res_10.00, fill = doublet_probability_median)) +
  geom_bar(stat = "identity") +
  labs(x= "Cluster", 
       y = "Median feature count",
       fill = "Doublet\nprobability") +
  my_theme() +
  theme(aspect.ratio = 1, legend.position = "inside", legend.position.inside = c(0.1,0.8))
dev.off()
