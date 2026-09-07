# Description -------------------------------------------------------------
# augur rank
# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)
iridescent <- khroma::color("iridescent")
# Load --------------------------------------------------------------------
full_results <- vroom::vroom(here::here("islet_cartography_scrna/data/cell_type_prioritization/files/augur_full_results.csv"))

# Preprocess --------------------------------------------------------------
# We first get the mean for each iteration to avoid pseudoreplication
input <- full_results |>
  tidyr::drop_na() |>
  dplyr::group_by(dataset, contrast, cell_type, idx) |>
  dplyr::summarise(
    mean_itx = mean(augur_score)) |>
  dplyr::ungroup() |> 
  dplyr::group_by(dataset, contrast, cell_type) |> 
  dplyr::summarise(
    mean = mean(mean_itx)) |> 
  dplyr::ungroup()


# Rank --------------------------------------------------------------------
scaled_rank <- input |> 
  dplyr::group_by(dataset, contrast) |> 
  dplyr::mutate(rank = base::rank(mean),
                rank_scaled = (rank - min(rank))/(max(rank) - min(rank))) |> 
  dplyr::ungroup() |> 
  dplyr::select(-mean, -rank) 

summary_rank <- scaled_rank |> 
  tidyr::pivot_wider(names_from = dataset, values_from = rank_scaled) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(mean_rank = mean(dplyr::c_across(cols = tidyselect::starts_with("ic_")), na.rm = T),
                median_rank = median(dplyr::c_across(cols = tidyselect::starts_with("ic_")), na.rm = T),
                sd_rank = sd(dplyr::c_across(cols = tidyselect::starts_with("ic_")), na.rm = T),
                sum_rank = sum(dplyr::c_across(cols = tidyselect::starts_with("ic_")), na.rm = T),
                n_not_na = sum(!is.na(dplyr::c_across(cols = tidyselect::starts_with("ic_"))))) |> 
  dplyr::select(!tidyselect::starts_with("ic_")) |> 
  dplyr::ungroup() |> 
  dplyr::filter(n_not_na > 4)


# Plot --------------------------------------------------------------------
pdf(here::here("islet_cartography_scrna/data/cell_type_prioritization/plot/mean_augur_score_nd_t2d.pdf"),
    width = 3,
    height = 2)
point_data <- scaled_rank |> 
  dplyr::filter(contrast == "nd_vs_t2d") |> 
  dplyr::filter(cell_type %in% c(summary_rank |> 
                                   dplyr::filter(contrast == "nd_vs_t2d") |> 
                                   dplyr::pull(cell_type) |> 
                                   unique()))

summary_rank |> 
  dplyr::filter(contrast == "nd_vs_t2d") |> 
  dplyr::arrange(mean_rank, dplyr::desc(sd_rank)) |> 
  dplyr::mutate(rank = dplyr::row_number()) |> 
  ggplot2::ggplot(ggplot2::aes(y = forcats::fct_reorder(cell_type, rank), x = mean_rank)) +
  ggplot2::geom_point(data = point_data, ggplot2::aes(y = cell_type, x = rank_scaled), 
                      color = "grey", size = 1) +
  ggplot2::geom_errorbar(aes(xmin = mean_rank - sd_rank, xmax = mean_rank + sd_rank), width = 0.2) +
  ggplot2::geom_point(ggplot2::aes(fill = n_not_na), size = 2, shape = 21) +
  ggplot2::geom_point(ggplot2::aes(x = median_rank), size = 2, shape= 21, 
                      color = "red", 
                      fill = NA) + 
  ggplot2::scale_fill_gradientn(colors = iridescent(5), limits=c(5, 11)) +
  ggplot2::labs(y = "Cell type",
                x = "Mean rank",
                fill = "Number of\n datasets") +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_title(x)) +
  my_theme() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank())
dev.off()


pdf(here::here("islet_cartography_scrna/data/cell_type_prioritization/plot/mean_augur_score_nd_pre.pdf"),
    width = 3,
    height = 2)
point_data <- scaled_rank |> 
  dplyr::filter(contrast == "nd_vs_pre") |> 
  dplyr::filter(cell_type %in% c(summary_rank |> 
                                   dplyr::filter(contrast == "nd_vs_pre") |> 
                                   dplyr::pull(cell_type) |> 
                                   unique()))

summary_rank |> 
  dplyr::filter(contrast == "nd_vs_pre") |> 
  dplyr::arrange(mean_rank, dplyr::desc(sd_rank)) |> 
  dplyr::mutate(rank = dplyr::row_number()) |> 
  ggplot2::ggplot(ggplot2::aes(y = forcats::fct_reorder(cell_type, rank), x = mean_rank)) +
  ggplot2::geom_point(data = point_data, ggplot2::aes(y = cell_type, x = rank_scaled), 
                      color = "grey", size = 1) +
  ggplot2::geom_errorbar(aes(xmin = mean_rank - sd_rank, xmax = mean_rank + sd_rank), width = 0.2) +
  ggplot2::geom_point(ggplot2::aes(fill = n_not_na), size = 2, shape = 21) +
  ggplot2::geom_point(ggplot2::aes(x = median_rank), size = 2, shape= 21, 
                      color = "red", 
                      fill = NA) + 
  ggplot2::scale_fill_gradientn(colors = iridescent(5), limits=c(5, 11)) +
  ggplot2::labs(y = "Cell type",
                x = "Mean rank",
                fill = "Number of\n datasets") +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_title(x)) +
  my_theme() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank())
dev.off()