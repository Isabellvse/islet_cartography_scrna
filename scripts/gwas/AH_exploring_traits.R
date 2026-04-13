# Description -------------------------------------------------------------
# Here I explore GWAS results

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load data ---------------------------------------------------------------
res <- vroom::vroom(here::here("islet_cartography_scrna/data/gwas/disease_scores/merged_scores.txt"))
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/meta.csv"))
t2d_nd <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/t2d_vs_nd.csv"))
nhood <- vroom::vroom(here::here("islet_cartography_scrna/data/milo/files/cells_in_neighborhood.csv"), altrep = TRUE)

gwas <- readxl::read_excel(here::here("islet_cartography_scrna/data/gwas/files/NIHMS1956161-supplement-Supp_Tables.xlsx"), 
                          sheet= "S1.Trait summary", skip = 1) |> 
  dplyr::rename_with(to_snake_case) |> 
  dplyr::select(trait, description, domain) |> 
  dplyr::distinct() |> 
  dplyr::mutate(category = dplyr::case_when(
    description %in% c(
      "Seen doctor (GP) for nerves, anxiety, tension or depression",
      "Fed-up feelings",
      "Guilty feelings",
      "Irritability",
      "Insomnia",
      "Balding Type 4",
      "Age at menarche",
      "Age at menopause",
      "Adult height"
    ) ~ "control",
    
    description %in% c(
      "Atrial fibrillation",
      "Asthma",
      "Autoimmune disease (Phecode + Self-reported)",
      "Breast cancer",
      "Coronary artery disease",
      "Cardiovascular disease",
      "Cholelithiasis",
      "Crohn's disease",
      "Diverticulitis",
      "Fibroblastic disorders",
      "Glaucoma (Phecode + Self-reported)",
      "Hypothyroidism",
      "Inflammatory bowel disease",
      "Lupus",
      "Blood clot in the leg"
    ) ~ "other_disease",
    
    description %in% c(
      "Type 2 diabetes",
      "Type 2 diabetes (adjusted by BMI)",
      "Body mass index",
      "Body fat percentage",
      "Body weight",
      "Glucose",
      "Hemoglobin A1c",
      "Hematocrit",
      "Hemoglobin",
      "Total protein",
      "Albumin",
      "Albumin/Globulin ratio",
      "Apolipoprotein A",
      "Apolipoprotein B",
      "High density lipoprotein cholesterol",
      "Low density lipoprotein cholesterol",
      "Total cholesterol",
      "Triglyceride",
      "C-reactive protein",
      "Calcium",
      "Vitamin D",
      "Uric acid",
      "Urea",
      "Total bilirubin",
      "Alkaline phosphatase",
      "Alanine aminotransferase",
      "Aspartate aminotransferase",
      "Gamma-glutamyl transferase",
      "Insulin-like growth factor 1",
      "Diastolic blood pressure",
      "Mean arterial pressure",
      "FEV1/FVC ratio",
      "Basophil count",
      "Eosinophil count",
      "Lymphocyte count",
      "White blood cell count",
      "Mean corpuscular hemoglobin",
      "Loss of Y"
    ) ~ "relevant",
    
    TRUE ~ NA_character_
  ))


# Define gradient colors
vik <- khroma::color("vik")
highcontrast <- color("high contrast")
# Preprocess --------------------------------------------------------------
df <- dplyr::full_join(res, meta)

## calculate median score - right now this is with norm score
df_median <- df |> 
  dplyr::group_by(manual_annotation, disease_harmonized) |> 
  dplyr::summarise(
    dplyr::across(tidyselect::ends_with("_norm"), list(median = median)), .groups = "drop"
  ) |> 
  tidyr::pivot_longer(cols = tidyselect::where(base::is.numeric), names_to = "trait", values_to = "median_score") |> 
  dplyr::mutate(trait = stringr::str_remove(trait, "_norm_median"), 
                cell_disease = base::paste0(manual_annotation, "_", disease_harmonized)) |> 
  # z-score
  dplyr::group_by(cell_disease) |> 
  dplyr::mutate(
    median_score_z = (median_score - mean(median_score, na.rm = TRUE)) / 
      sd(median_score, na.rm = TRUE)
  ) |> 
  dplyr::ungroup() |> 
  dplyr::select(-median_score)


# Plot as heatmap ---------------------------------------------------------
## Matrix
mat <- df_median |>
  tidyr::pivot_wider(
    names_from  = trait,
    values_from = median_score_z
  ) |>
  tibble::column_to_rownames("cell_disease") |>
  dplyr::select(tidyselect::where(is.numeric)) |> 
  as.matrix()

## Hierarchical Clustering
row_hc <- stats::hclust(stats::dist(mat, method = "euclidean"), method = "complete")
col_hc <- stats::hclust(stats::dist(base::t(mat), method = "euclidean"), method = "complete")

## Set row and col order
row_order <- rownames(mat)[row_hc$order]
col_order <- colnames(mat)[col_hc$order]

## Add trait description
cell_disease_order <- df_median_clust$cell_disease |> unique()

df_median_clust <- df_median |>
  dplyr::left_join(y = gwas, by = "trait") |> 
  dplyr::mutate(trait = base::factor(trait, levels = col_order),
                cell_disease = base::factor(cell_disease, levels = cell_disease_order))

## Range of data
rng <- range(df_median_clust$median_score_z, na.rm = TRUE)

### Plot heatmap with clipped values
# everything above 1 or below -1 is set to 1 or -1 respectively


# p_score <- df_median_clust |>
#   ggplot2::ggplot(
#     ggplot2::aes(y = trait, x = cell_disease, fill = median_score_z)
#   ) +
#   ggplot2::geom_tile(color = "darkgrey") +
#   ggplot2::scale_fill_gradientn(
#     colors = vik(256),
#     limits = c(-2, 2),
#     values = scales::rescale(c(-2, 0, 2)),
#     oob = scales::squish
#   ) +
#   my_theme() +
#   ggplot2::theme(
#     axis.text.y = ggplot2::element_blank(),
#     axis.ticks.y = ggplot2::element_blank(),
#     axis.title.y = ggplot2::element_blank(),
#     axis.line.y = ggplot2::element_blank(),
#     axis.text.x = ggplot2::element_blank(),
#     axis.ticks.x = ggplot2::element_blank(),
#     axis.title = ggplot2::element_blank()
#   )


p_score <- df_median_clust |>
  ggplot2::ggplot(
    ggplot2::aes(y = trait, x = cell_disease, fill = median_score_z)
  ) +
  ggplot2::geom_tile(color = "darkgrey") +
  ggplot2::scale_fill_gradientn(
    colors = vik(256),
    values = scales::rescale(c(rng[1], 0, rng[2])),
    limits = rng
  ) +
  my_theme() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank()
  )

p_domain <- df_median_clust |>
  ggplot2::ggplot(
    ggplot2::aes(y = trait, x = "category", fill = category)
  ) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_manual(values = c(highcontrast(3), "grey"))+
  my_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank())

p_celltype <- df_median_clust |>
  ggplot2::ggplot(
    ggplot2::aes(y = "1", x = cell_disease, fill = manual_annotation)
  ) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_manual(values = manual_anno_colors)+
  my_theme() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))


p_disease <- df_median_clust |>
  ggplot2::ggplot(
    ggplot2::aes(y = "1", x = cell_disease, fill = disease_harmonized)
  ) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_manual(values = disease_color)+
  my_theme() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank())


# Top row: p_domain + p_score
top_row <- p_domain + p_score + plot_layout(widths = c(0.1, 5))

# Bottom row: spacer under p_domain + p_disease / p_celltype under p_score
bottom_row <- plot_spacer() + (p_disease / p_celltype) + plot_layout(widths = c(0.1, 5), heights = c(1,1))

# Combine top and bottom rows
p_combined <- top_row / bottom_row + plot_layout(heights = c(5, 2), guides = "collect")

p_combined


# t2d hoods ---------------------------------------------------------------
beta_sig <- t2d_nd |> dplyr::filter(FDR <= 0.10 & nhood_annotation == "beta")
nhood_sig <- nhood |> dplyr::filter(cell %in% unique(beta_sig$index_cell)) |> 
  dplyr::distinct()


test_2 <- test |> dplyr::filter(padj <= 0.05) 


df_sig <- df |> dplyr::filter(barcode %in% unique(beta_sig$index_cell)) |> 
  ggplot2::ggplot(aes(y = FDR, x=logFC)) 
  
test<-vroom::vroom(here::here("islet_cartography_scrna/data/annotate/dg_onevsother/deg_wald_manual_annotation_endmt_early_vs_allother.csv"))
test_2 <- test |> dplyr::filter(padj <= 0.05 & pct_expr_c1 > 60) |> top_n(n = 20, wt = log2FoldChange)
