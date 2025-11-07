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


# Evaluate doublet threshold per cluster ----------------------------------
# threshold 0.4
c04 <- adata_sum  |>  
  dplyr::mutate(filter = dplyr::case_when(doublet_probability_median >= 0.4 ~ 1,
                                          .default = 0)) |> 
  dplyr::filter(filter >= 1) |> 
  dplyr::pull(leiden_res_10.00)

# threshold 0.5
c05 <- adata_sum  |>  
  dplyr::mutate(filter = dplyr::case_when(doublet_probability_median >= 0.5 ~ 1,
                                          .default = 0)) |> 
  dplyr::filter(filter >= 1) |> 
  dplyr::pull(leiden_res_10.00)

# Find difference between 0.4 and 0.5
diff_clusters <- setdiff(c04, c05)


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
"#d3d3d3", "black"

## Plot histograms ---------------------------------------------------------
# Median module score
adata_comb |> 
  dplyr::pull(module_difference) |> 
  plot_hist_qc(.y = module_difference, title = "Median module score difference")


# Evaluate new doublet threshold ------------------------------------------
# Evaluate doublet threshold per cluster ----------------------------------
# threshold 0.4
c042 <- adata_comb  |>  
  dplyr::mutate(db_filter = dplyr::case_when(doublet_probability_median >= 0.4 ~ 1,
                                          .default = 0),
                mo_filter = dplyr::case_when(module_difference <= 2 ~ 1,
                                          .default = 0),
                filter = db_filter + mo_filter) |> 
  dplyr::filter(filter >= 2) |> 
  dplyr::pull(leiden_res_10.00)

# threshold 0.5
c052 <- adata_comb  |>  
  dplyr::mutate(db_filter = dplyr::case_when(doublet_probability_median >= 0.5 ~ 1,
                                             .default = 0),
                mo_filter = dplyr::case_when(module_difference <= 2 ~ 1,
                                             .default = 0),
                filter = db_filter + mo_filter) |> 
  dplyr::filter(filter >= 2) |> 
  dplyr::pull(leiden_res_10.00)

# Find difference between 0.4 and 0.5
diff_clusters2 <- setdiff(c042, c052)

diff_clusters
diff_clusters2


# Diff genes --------------------------------------------------------------
# Add identifier for split samples by 
test <- adata_obs |> 
  dplyr::select(ic_id_dataset_donor, ic_id_sample, platform) |> 
  dplyr::mutate(ic_id_dataset_donor_sample = 
                  dplyr::case_when(platform %in% c("droplet", "plate_barcode") ~ ic_id_sample,
                                   platform %in% c("plate") ~ ic_id_dataset_donor)) |> 
  dplyr::distinct() 


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

adata_sum  |> 
  dplyr::filter(doublet_probability_median >= 0.5) |> 
  dplyr::pull(leiden_res_10.00)

adata_sum  |>  
  dplyr::filter(leiden_res_10.00 %in% filter_cluster) |> 
  mutate(leiden_res_10.00 = forcats::fct_reorder(as.character(leiden_res_10.00), dplyr::desc(n_cells))) |> 
  ggplot2::ggplot(aes(x = n_cells, y = leiden_res_10.00)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_text(aes(label = n_cells), hjust = 0, size = 2) +
  ggplot2::labs(title = "Cluster size",
                y = "Leiden cluster res 10.00",
                x = "Number of cells ") +
  my_theme()


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


## Define doublets ----
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

plot_list[[1]] + ggplot2::geom_vline(xintercept = 2000, color = "red") +
plot_list[[6]] + ggplot2::geom_vline(xintercept = 0.4, color = "red") +
plot_list[[7]] + ggplot2::geom_vline(xintercept = 0.05, color = "red") +
plot_list[[8]] + ggplot2::geom_vline(xintercept = 0.5, color = "red")




# A closer look on clusters in the middle ---------------------------------
test <- 

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
# Per cell analysis -------------------------------------------------------
# Dependencies
library(lme4)
library(lmerTest)

# Helper function ----------------------------------------------------------
analyze_metric <- function(data, value_col, label, round_digits = 0) {
  # Prepare data
  agg_df <- data %>%
    filter(!is.na(n_umi)) %>%
    group_by(ic_id_donor_overall, ic_id_sample, doublet) %>%
    summarise(median_value = median(.data[[value_col]]), .groups = "drop")
  
  # Fit mixed model
  model <- lmer(median_value ~ doublet + (1 | ic_id_donor_overall), data = agg_df)
  fe <- summary(model)$coefficients["doubletsinglet", ]
  
  estimate <- round(fe["Estimate"], round_digits)
  p_value <- signif(fe["Pr(>|t|)"], 2)
  t_value <- fe["t value"]
  df <- fe["df"]
  
  # Effect size (Rosenthal 1991, p 19) 
  # https://peterstatistics.com/CrashCourse/3-TwoVarUnpair/BinOrd/BinOrd-2b-EffectSize.html
  r <- sqrt(t_value^2 / (t_value^2 + df))
  
  # Annotation text
  annotation_text <- paste0(
    "Linear mixed model:\n",
    "Delta ", label, " ~= ", estimate,
    "\np = ", p_value,
    "\nEffect size (Rosenthal) = ", round(r, 2)
  )
  labs(
    x = "Doublet classification\n(>= 90% probability)"
  )
  
  # y-position
  y_pos <- max(agg_df$median_value, na.rm = TRUE) * 0.90
  
  # Plot
  p <- agg_df %>%
    mutate(doublet = factor(doublet, levels = c("singlet", "doublet"))) %>%
    ggplot(aes(x = doublet, y = median_value)) +
    geom_boxplot(width = 0.2, fill = "transparent", outlier.size = 0.0001, outlier.color = "grey") +
    geom_text(
      data = data.frame(x = 0.5, y = y_pos, label = annotation_text),
      aes(x = x, y = y, label = label),
      size = 1,
      hjust = 0  # left-align
    ) +
    labs(
      title = label,
      x = paste0("Doublet classification\n(>= 90% probability)"),
      y = paste0("Median ", label, "\nper sample")
    ) +
    my_theme() +
    theme(aspect.ratio = 1)
  
  return(p)
}


## Prepare metrics ---------------------------------------------------------
max_score_col <- adata_obs %>%
  dplyr::select(barcode, tidyselect::matches("^azimuth.*_score")) %>%
  tidyr::pivot_longer(cols = -barcode, names_to = "score_col", values_to = "value") %>%
  dplyr::group_by(barcode) %>%
  dplyr::slice_max(order_by = value, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(barcode, score_col)

# Join max score info back to adata_sum
adata_sum_max <- adata_obs %>%
  dplyr::left_join(max_score_col, by = "barcode")

# Calculate absolute differences from the max score
adata_sum_dist <- adata_sum_max %>%
  dplyr::rowwise() %>%
  dplyr::mutate(across(
    tidyselect::matches("^azimuth.*_score"),
    ~ abs(. - get(score_col)),
    .names = "{.col}_dist_from_max"
  )) %>%
  dplyr::ungroup()

test <- adata_sum_dist |>
  dplyr::select(barcode, tidyselect::contains("_dist_from_max"), score_col) |> 
  tidyr::pivot_longer(c(-barcode, -score_col, -score_col)) |> 
  dplyr::group_by(barcode) |> 
  dplyr::mutate(remove = case_when(value == 0 & score_col == stringr::str_remove(name, "_dist_from_max") ~ TRUE,
                                   .default = FALSE)) |> 
  dplyr::rename(module_name = name) |> 
  dplyr::filter(!remove == TRUE) |> 
  dplyr::select(-remove) |> 
  dplyr::slice_min(order_by = value, n= 1) |> 
  dplyr::ungroup() |> 
  dplyr::right_join(adata_obs) |> 
  dplyr::mutate(doublet = dplyr::case_when(doublet_probability >= 0.9 ~ "doublet",
                                           doublet_probability < 0.9 ~ "singlet",
                                           .default = "singlet"))

## Run analysis for each metric --------------------------------------------
p_umi <- analyze_metric(test, "n_umi", "UMI", 0)
p_diff <- analyze_metric(test, "value", "Module score difference", 2)
p_feature <- analyze_metric(test, "n_feature", "Feature", 0)

## Save combined plot ------------------------------------------------------
pdf(
  file = here::here("islet_cartography_scrna/data/integrate/second_pass/plot/doublet_boxplot.pdf"),
  height = 2,
  width = 5
)
p_diff + p_umi + p_feature
dev.off()


# Per cluster analysis ----------------------------------------------------

# Preprocess --------------------------------------------------------------
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


# Distances in gene set scores --------------------------------------------
# Identify which leiden cluster that has the max score
max_score_col <- adata_sum %>%
  dplyr::select(leiden_res_10.00, tidyselect::matches("^azimuth.*_score_median$")) %>%
  tidyr::pivot_longer(cols = -leiden_res_10.00, names_to = "score_col", values_to = "value") %>%
  dplyr::group_by(leiden_res_10.00) %>%
  dplyr::slice_max(order_by = value, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(leiden_res_10.00, score_col)

# Join max score info back to adata_sum
adata_sum_max <- adata_sum %>%
  dplyr::left_join(max_score_col, by = "leiden_res_10.00")

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
  dplyr::select(leiden_res_10.00, tidyselect::contains("_dist_from_max"), score_col) |> 
  tidyr::pivot_longer(c(-leiden_res_10.00, -score_col, -score_col), ) |> 
  dplyr::group_by(leiden_res_10.00) |> 
  dplyr::mutate(remove = case_when(value == 0 & score_col == stringr::str_remove(name, "_dist_from_max") ~ TRUE,
                                   stringr::str_remove(name, "_dist_from_max") == "azimuth_cycling_score_median" ~ TRUE,
                                   .default = FALSE)) |> 
  dplyr::filter(!remove == TRUE) |> 
  dplyr::select(-remove) |> 
  dplyr::slice_min(order_by = value, n= 1) |> 
  dplyr::ungroup() |> 
  dplyr::right_join(adata_obs |> dplyr::select(leiden_res_10.00, doublet_probability)) |> 
  dplyr::mutate(doublet = dplyr::case_when(doublet_probability >= 0.9 ~ "doublet",
                                           doublet_probability < 0.9 ~ "singlet",
                                           .default = "singlet")) |> 
  dplyr::group_by(leiden_res_10.00) |> 
  dplyr::mutate(perc_doublet = sum(doublet == "doublet", na.rm = TRUE) / dplyr::n() * 100) |> 
  dplyr::ungroup()

# Extract cluster based metrics
test_2 <- test |> 
  dplyr::select(-doublet_probability, -doublet) |> dplyr::distinct()

# Plot perc of doublets vs minimum distance in gene module scores
plot(test_2$perc_doublet, test_2$value)
abline(h=0.50, v =10, col = "red")

# Try to set a threshold, 
# those clusters with a distance less than 0.8 and a doublet % higher than 20
# are marked as doublet clusters
doublet_clusters <- test_2 |> 
  dplyr::filter(value <= 0.50 & perc_doublet >= 10) |> 
  pull(leiden_res_10.00)

test_2 <- adata_sum_dist |>
  dplyr::select(leiden_res_10.00, tidyselect::contains("_dist_from_max"), score_col) |> 
  tidyr::pivot_longer(c(-leiden_res_10.00, -score_col, -score_col), ) |> 
  dplyr::group_by(leiden_res_10.00) |> 
  dplyr::mutate(remove = case_when(value == 0 & score_col == stringr::str_remove(name, "_dist_from_max") ~ TRUE,
                                   stringr::str_remove(name, "_dist_from_max") == "azimuth_cycling_score_median" ~ TRUE,
                                   .default = FALSE)) |> 
  dplyr::filter(!remove == TRUE) |> 
  dplyr::select(-remove) |> 
  dplyr::slice_min(order_by = value, n= 1) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(doublet = dplyr::case_when(leiden_res_10.00 %in% doublet_clusters ~ "doublet",
                                           .default = "singlet")) |> 
  dplyr::right_join(adata_sum_dist)


# Helper function ----------------------------------------------------------
analyze_metric <- function(data, value_col, label, round_digits = 0) {
  # Prepare data
  agg_df <- data |> 
    dplyr::rename(median_value = value_col)

  # Plot
  p <- agg_df %>%
    mutate(doublet = factor(doublet, levels = c("singlet", "doublet"))) %>%
    ggplot(aes(x = doublet, y = median_value)) +
    geom_boxplot(width = 0.2, fill = "transparent", outlier.size = 0.0001, outlier.color = "grey") +

    labs(
      title = label,
      x = paste0("Doublet classification"),
      y = paste0("Median ", label, "\nper leiden cluster")
    ) +
    my_theme() +
    theme(aspect.ratio = 1)
  
  return(p)
}


## Run analysis for each metric --------------------------------------------
p_umi <- analyze_metric(test_2, "n_umi_median", "UMI", 0)
p_diff <- analyze_metric(test_2, "value", "Module score difference", 2)
p_feature <- analyze_metric(test_2, "n_feature_median", "Feature", 0)

pdf(
  file = here::here("islet_cartography_scrna/data/integrate/second_pass/plot/doublet_boxplot_per_cluster.pdf"),
  height = 2,
  width = 5
)
p_diff + p_umi + p_feature
dev.off()

# Other quality control per cluster ---------------------------------------

plot_hist_qc <- function(.x, .y, title) {
  ggplot2::ggplot(data = tibble::tibble(value = .x), ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 50, fill = "black") +
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
  file = here::here("islet_cartography_scrna/data/integrate/second_pass/plot/qc_plots_per_cluster_leiden_5.pdf"),
  height = 2,
  width = 2
)
adata_sum |> 
  dplyr::select(tidyselect::contains("_median"), -tidyselect::contains("azimuth")) |>
  purrr::iwalk(~ print(plot_hist_qc(.x, .y, title = .y)))
dev.off()




pdf(width = 12, height = 2)
wrap_plots(plot_list, nrow = 1) + plot_annotation(title = "Resolution = 10.00")
dev.off()

pdf(width = 12, height = 2)
wrap_plots(plot_list, nrow = 1) + plot_annotation(title = "Resolution = 5.50")
dev.off()

pdf(width = 12, height = 2)
wrap_plots(plot_list_2, nrow = 1) + plot_annotation(title = "Resolution = 10.00")
dev.off()


plot_list <-adata_sum %>% 
  dplyr::select(tidyselect::contains("_median"), -tidyselect::contains("azimuth")) |>
  purrr::imap(~ plot_hist_qc(.x, .y, title = .y))


plot_list_2 <- plot_list %>% 
  purrr::imap(function(.x, .y){
    vec <- adata_sum %>% dplyr::filter(contrast_fraction_median < 0.7) %>% dplyr::pull(.y) 
    output <- .x + ggplot2::geom_vline(xintercept = vec, color = "red") 
    return(output)
  })

clust <- adata_sum %>% dplyr::filter(contrast_fraction_median < 0.7) %>% dplyr::pull("leiden_res_10.00")

adata_obs %>% filter(leiden_res_10.00 %in% clust) %>% group_by(cell_nuclei) %>% tally()
adata_obs %>% filter(!leiden_res_10.00 %in% clust) %>% group_by(cell_nuclei) %>% tally()

umi_vec <- adata_sum %>% dplyr::filter(contrast_fraction_median < 0.7) %>% dplyr::pull(n_umi_median) 

plot_list[[1]] & ggplot2::geom_vline(xintercept = umi_vec, color = "red") 
