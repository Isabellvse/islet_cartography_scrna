# Description -------------------------------------------------------------
# Cell type proportion annotations using binomial mixed models
# Cell type proportions were modelled using a binomial generalised linear mixed model (GLMM) 
# implemented in glmmTMB. For each cell type, the number of cells of that type per sample was 
# modelled against the total number of cells in that sample using cbind(n_cells, total - n_cells), 
# which specifies the number of successes (cells of the target type) and failures (all other cells) respectively. 

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)
library(glmmTMB)
library(emmeans)

create_directories(list(here::here("islet_cartography_scrna/data/cell_type_proportions"),
                        here::here("islet_cartography_scrna/data/cell_type_proportions/plot"),
                        here::here("islet_cartography_scrna/data/cell_type_proportions/files")))

# Load --------------------------------------------------------------------
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/obs.csv"))
duplicated_donors <- vroom::vroom(here::here("islet_cartography_scrna/data/duplicated_donors.csv"))

# Preprocess - filter donors ----------------------------------------------
# Keep single-dataset donors + 10x versions of multi-dataset donors
meta_filtered <- meta |> 
  dplyr::anti_join(y = duplicated_donors)

stopifnot(
  meta_filtered |> 
    dplyr::distinct(ic_id_donor_overall, ic_id_dataset) |> 
    dplyr::add_count(ic_id_donor_overall) |> 
    dplyr::pull(n) |> 
    max() == 1
)

print(glue::glue("Removed {nrow(meta) - nrow(meta_filtered)} observations"))

# preprocess data for model -----------------------------------------------
# Number of cell per donor per study
cell_counts <- meta_filtered |> 
  dplyr::group_by(ic_id_dataset, ic_id_donor_overall, cell_type_broad) |>
  dplyr::summarise(n_cells = dplyr::n(), .groups = 'drop') |>
  dplyr::group_by(ic_id_dataset, ic_id_donor_overall) |>  
  dplyr::mutate(total = sum(n_cells)) |>
  dplyr::ungroup()

# donor level meta data
meta_counts <- meta_filtered |>
  dplyr::select(ic_id_donor_overall, ic_id_dataset,
                disease_hba1c, bmi, sex_predicted, age_years) |>
  dplyr::distinct() |>
  dplyr::left_join(
    cell_counts,
    by = c("ic_id_dataset", "ic_id_donor_overall")) |>
  # Remove donors with missing values
  dplyr::filter(!is.na(bmi), !is.na(age_years), !is.na(sex_predicted)) |>
  dplyr::mutate(
    disease_hba1c = factor(disease_hba1c, levels = c("nd", "pre", "t2d")),
    sex_predicted      = factor(sex_predicted, levels = c("m", "f")),
    ic_id_dataset      = factor(ic_id_dataset),
    age_z              = as.numeric(scale(age_years)),   
    bmi_z              = as.numeric(scale(bmi))
  )

# glmm --------------------------------------------------------------------
fits <- meta_counts |> 
  base::split(~cell_type_broad) |> 
  purrr::map(\(df) {glmmTMB::glmmTMB(
      as.formula("cbind(n_cells, total - n_cells) ~ disease_hba1c + age_z + sex_predicted + bmi_z + (1| ic_id_dataset)"), 
      family = betabinomial, data = df)})

coef_results <- fits |> 
  purrr::map(\(fit){
    fit |> broom.mixed::tidy(effects = "fixed")
  })

# estimated means ---------------------------------------------------------
emm_mean <- fits |> 
  purrr::map(\(fit) {emmeans::emmeans(fit, ~ disease_hba1c, 
                     type = "response", 
                     weights = "proportional") |> 
      summary() |> 
      as.data.frame()
  }) |> 
  purrr::list_rbind(names_to = "cell_type_broad")

emm_res <- fits |> 
  purrr::map(\(fit){
    emmeans::emmeans(fit, pairwise ~ disease_hba1c, 
                     adjust = "none",
                     type = "response", 
                     weights = "proportional")$contrasts  |> 
      as.data.frame()
  }) |> 
  purrr::list_rbind(names_to = "cell_type_broad") |> 
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr"), 
                star = dplyr::if_else(fdr <= 0.05, "*", ""))

# Plot --------------------------------------------------------------------
# get max mean value
emm_max <- emm_mean |> 
  dplyr::group_by(cell_type_broad) |> 
  dplyr::summarise(ymax = max(asymp.UCL), .groups = "drop")

# get significant values
emm_sig <- emm_res |> 
  dplyr::filter(star == "*") |>          # only plot significant ones
  dplyr::left_join(emm_max) |> 
  dplyr::group_by(cell_type_broad) |> 
  dplyr::mutate(y.position = ymax * (1 + 0.12 * dplyr::row_number())) |> 
  tidyr::separate(contrast, into = c("group1", "group2"), sep = " / ") |>
  dplyr::ungroup()

pdf(here::here(glue::glue("islet_cartography_scrna/data/cell_type_proportions/plot/cell_type_broad_proportions.pdf")),
    height = 1.5, 
    width = 7)
emm_mean |>
  ggplot2::ggplot(ggplot2::aes(x = disease_hba1c, y = prob, fill = cell_type_broad)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.3) +
  ggpubr::stat_pvalue_manual(
    emm_sig,
    label = "star",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = NA
  ) +
  ggplot2::facet_wrap(~cell_type_broad, scales = "free_y", nrow = 1, labeller = labeller(
    cell_type_broad = function(x)
      ggplot2::label_wrap_gen(width = 10)(stringr::str_to_title(gsub("_", " ", x))
      ))) +
  ggplot2::scale_x_discrete(
    labels = function(x)(
      stringr::str_to_upper(x))
  ) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ggplot2::labs(x = "Disease status",
                y = "Estimated proportion") +
  ggplot2::scale_fill_manual(values = cell_type_broad_colors) +
  my_theme() +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


# Subtypes ----------------------------------------------------------------
subtype_counts <- meta_filtered |>
  dplyr::group_by(
    ic_id_dataset,
    ic_id_donor_overall,
    cell_type_broad,
    cell_type
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::group_by(
    ic_id_dataset,
    ic_id_donor_overall,
    cell_type_broad
  ) |>
  dplyr::mutate(
    total_parent = sum(n_cells)
  ) |>
  dplyr::ungroup()

subtype_counts <- meta_filtered |>
  dplyr::select(ic_id_donor_overall, ic_id_dataset,
                disease_hba1c, bmi, sex_predicted, age_years) |>
  dplyr::distinct() |>
  dplyr::left_join(
    subtype_counts,
    by = c("ic_id_dataset", "ic_id_donor_overall")) |>
  # Remove donors with missing values
  dplyr::filter(!is.na(bmi), !is.na(age_years), !is.na(sex_predicted)) |>
  dplyr::mutate(
    disease_hba1c = factor(disease_hba1c, levels = c("nd", "pre", "t2d")),
    sex_predicted      = factor(sex_predicted, levels = c("m", "f")),
    ic_id_dataset      = factor(ic_id_dataset),
    age_z              = as.numeric(scale(age_years)),   
    bmi_z              = as.numeric(scale(bmi))
  ) %>% 
  dplyr::filter(cell_type_broad %in% c("stellate", "acinar", "ductal"))

fits_subtype <- subtype_counts |>
  base::split(~cell_type) |>
  purrr::map(\(df) {glmmTMB(
    cbind(n_cells, total_parent - n_cells) ~
      disease_hba1c +
      age_z +
      sex_predicted +
      bmi_z +
      (1| ic_id_dataset),
    family = betabinomial(),
    data = df)})

# estimated means ---------------------------------------------------------
emm_mean <- fits_subtype  |> 
  purrr::map(\(fit){
    emmeans::emmeans(fit, ~ disease_hba1c, 
                     type = "response",
                     weights = "proportional") |> 
      summary() |> 
      as.data.frame()
  }) |> 
  purrr::list_rbind(names_to = "cell_type")

emm_res <- fits_subtype  |> 
  purrr::map(\(fit){
    emmeans::emmeans(fit, pairwise ~ disease_hba1c, 
                     adjust = "none",
                     type = "response",
                     weights = "proportional")$contrasts  |> 
      as.data.frame()
  }) |> 
  purrr::list_rbind(names_to = "cell_type") |> 
  dplyr::mutate(fdr = stats::p.adjust(p.value, method = "fdr"), 
                star = dplyr::if_else(fdr <= 0.05, "*", ""))
# Plot --------------------------------------------------------------------
# get max mean value
emm_max <- emm_mean |> 
  dplyr::group_by(cell_type) |> 
  dplyr::summarise(ymax = max(asymp.UCL), .groups = "drop")

# get significant values
emm_sig <- emm_res |> 
  dplyr::filter(star == "*") |>          # only plot significant ones
  dplyr::left_join(emm_max) |> 
  dplyr::group_by(cell_type) |> 
  dplyr::mutate(y.position = ymax * (1 + 0.12 * dplyr::row_number())) |> 
  tidyr::separate(contrast, into = c("group1", "group2"), sep = " / ") |>
  dplyr::ungroup()

pdf(here::here(glue::glue("islet_cartography_scrna/data/cell_type_proportions/plot/cell_type_subtype_proportions.pdf")),
    height = 1.5, 
    width = 2.7)
emm_mean |>
  ggplot2::ggplot(ggplot2::aes(x = disease_hba1c, y = prob, fill = cell_type)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.3) +
  # ggpubr::stat_pvalue_manual(
  #   emm_sig,
  #   label = "star",
  #   xmin = "group1", xmax = "group2",
  #   y.position = "y.position",
  #   tip.length = NA
  # ) +
  ggplot2::facet_wrap(~cell_type, scales = "free_y", nrow = 1, labeller = labeller(
    cell_type = function(x)
      ggplot2::label_wrap_gen(width = 10)(stringr::str_to_title(gsub("_", " ", x))
      ))) +
  ggplot2::scale_x_discrete(
    labels = function(x)(
      stringr::str_to_upper(x))
  ) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ggplot2::labs(x = "Disease status",
                y = "Estimated cell-type proportion") +
  ggplot2::scale_fill_manual(values = cell_type_colors) +
  my_theme() +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

