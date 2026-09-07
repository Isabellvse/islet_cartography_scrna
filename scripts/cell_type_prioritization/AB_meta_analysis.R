# Description -------------------------------------------------------------
# augur rank
# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
full_results <- vroom::vroom(here::here("islet_cartography_scrna/data/cell_type_prioritization/files/augur_full_results.csv"))

# Preprocess --------------------------------------------------------------
# We first get the mean for each iteration to avoid
# pseudoreplication
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
test <- input %>% 
  dplyr::group_by(dataset, contrast) |> 
  dplyr::mutate(rank = base::rank(mean),
                rank_scaled = (rank - min(rank))/(max(rank) - min(rank))) |> 
  dplyr::ungroup() |> 
  dplyr::select(-mean, -rank) |> 
  tidyr::pivot_wider(names_from = dataset, values_from = rank_scaled) |> 
  dplyr::rowwise() |> 
  dplyr::mutate(mean_rank = mean(dplyr::c_across(cols = tidyselect::starts_with("ic_")), na.rm = T),
                median_rank = median(dplyr::c_across(cols = tidyselect::starts_with("ic_")), na.rm = T),
                sd_rank = sd(dplyr::c_across(cols = tidyselect::starts_with("ic_")), na.rm = T),
                sum_rank = sum(dplyr::c_across(cols = tidyselect::starts_with("ic_")), na.rm = T),
                n_not_na = sum(!is.na(dplyr::c_across(cols = tidyselect::starts_with("ic_"))))) |> 
  dplyr::select(!tidyselect::starts_with("ic_")) |> 
  dplyr::filter(n_not_na > 4)
