# Description -------------------------------------------------------------
# Plotting meta data from donors left from second around of QC (per cluster)
# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)


# Load --------------------------------------------------------------------
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/post_qc_metadata.csv"))


# preprocess --------------------------------------------------------------

### Meta data per donor ----
meta_donor <- meta |> 
  dplyr::select(name, ic_id_donor_overall, ic_id_platform_adjusted_sample, ic_id_study, ic_id_dataset, library_prep, disease_harmonized,
                sex_predicted, age_years, bmi, hba_1_c_percent, ethnicity_broad_harmonized, diabetes_medication_harmonized, 
                cause_of_death_broad_harmonized, cell_nuclei) |> 
  dplyr::distinct()

# Statistical tests -------------------------------------------------------
## blood glucose ----
df <- meta_donor |> 
  dplyr::select(ic_id_donor_overall, hba_1_c_percent, disease_harmonized) |> 
  dplyr::distinct() |> 
  tidyr::drop_na()

# Fit ANOVA model
aov_model <- aov(hba_1_c_percent ~ disease_harmonized, data = df)

# Tukey post-hoc test
rstatix::tukey_hsd(aov_model)

# Effect size
df |> 
  rstatix::cohens_d(hba_1_c_percent ~ disease_harmonized)


# heatmap -----------------------------------------------------------------
cols22 <- c(
  "#1F51FF", "#d95f02", "#CBC3E3", "#e7298a", "#66a61e",
  "#e6ab02", "#800020", "#988558",
  "#00FFFF", "#b2df8a", "#fb9a99", "#FFEA00",
  "#FF00FF", "#5D3FD3", "#0FFF50", "#b15928",
  "#8dd3c7", "#bebada", "#80b1d3", "#fccde5",
  "#b3de69", "#FF5F1F"
)

cols25 <- c(
  "#1F51FF", "#d95f02", "#CBC3E3", "#e7298a", "#66a61e",
  "#e6ab02", "#800020", "#988558",
  "#00FFFF", "#b2df8a", "#fb9a99", "#FFEA00",
  "#FF00FF", "#5D3FD3", "#0FFF50", "#b15928",
  "#8dd3c7", "#bebada", "#80b1d3", "#fccde5",
  "#b3de69", "#FF5F1F", "#FFA82E", "#DFFF00", "#FAD5A5"
)

cols5 <- c("#1F51FF", "#FF00FF", "#0FFF50", "#FF5F1F", "#FFEA00")

# Donor order
donor_order <- meta_donor |>
  dplyr::arrange(
    ic_id_study,
    disease_harmonized,
    hba_1_c_percent,
    ic_id_dataset,
  ) |>
  dplyr::distinct(ic_id_donor_overall, .keep_all = TRUE) |>
  dplyr::pull(ic_id_donor_overall)

tile_row <- function(data, fill_var, y_label) {
  data |>
    distinct(ic_id_donor_overall, {{ fill_var }}) |>
    mutate(ic_id_donor_overall = factor(ic_id_donor_overall, levels = donor_order)) |>
    ggplot(aes(x = ic_id_donor_overall, y = 1, fill = {{ fill_var }})) +
    geom_tile(color = NA) +
    labs(y = y_label, fill = NULL) +
    my_theme() +
    theme(
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "none"
    )
}

study   <- tile_row(meta_donor, ic_id_study,   "Study") + scale_fill_manual(values = cols22)
disease <- tile_row(meta_donor, disease_harmonized, "Disease") + ggplot2::scale_fill_manual(values = disease_color)
sex <- tile_row(meta_donor, sex_predicted, "Sex") + scale_fill_manual(values = c("#C2563A", "#3F7F93"))
eth <- tile_row(meta_donor, ethnicity_broad_harmonized, "Ethnicity") + scale_fill_manual(values = cols5)
age  <- tile_row(meta_donor, age_years, "Age (years)") + ggplot2::scale_fill_gradient(low = "white", high = "#5D3FD3")
bmi  <- tile_row(meta_donor, bmi, "BMI") + ggplot2::scale_fill_gradient(low = "white", high = "#FF5F1F")
medication <- tile_row(meta_donor, diabetes_medication_harmonized, "Diabetes medicaiton") + scale_fill_manual(values = c("white", "red"))

dataset <-  meta_donor |>
  count(ic_id_donor_overall, ic_id_dataset) |>
  mutate(ic_id_donor_overall = factor(ic_id_donor_overall, levels = donor_order)) |>
  ggplot(aes(x = ic_id_donor_overall, y = n, fill = ic_id_dataset)) +
  geom_col(position = "fill", width = 1, color = NA, linewidth = 0.3) +
  labs(y = "Dataset", fill = NULL) +
  scale_fill_manual(values = cols25) +
  my_theme() +
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  )

platform <- meta_donor |>
  count(ic_id_donor_overall, library_prep) |>
  mutate(ic_id_donor_overall = factor(ic_id_donor_overall, levels = donor_order)) |>
  ggplot(aes(x = ic_id_donor_overall, y = n, fill = library_prep)) +
  geom_col(position = "fill", width = 1, color = NA, linewidth = 0.3) +
  labs(y = "Platform", fill = NULL) +
  scale_fill_manual(values = cols25) +
  my_theme() +
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  )

hba <- meta_donor |>
  distinct(ic_id_donor_overall, hba_1_c_percent) |>
  mutate(ic_id_donor_overall = factor(ic_id_donor_overall, levels = donor_order)) |>
  ggplot(aes(x = ic_id_donor_overall, y = 1, fill = hba_1_c_percent)) +
  geom_tile(color = NA) +
  scale_fill_gradientn(
    colours = c("forestgreen", "white", "gold", "#ae0000"),
    values = scales::rescale(c(4.5, 5.7, 6.5, 13.1)),
    breaks = c(4.5, 5.7, 6.5, 13.1),
    name = NULL
  ) +
  labs(y = "HbA1c %") +
  my_theme() +
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
    legend.position = "none"
  )


# get legends -------------------------------------------------------------
library(cowplot)

get_leg <- function(p) {
  cowplot::get_legend(
    p +
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(nrow = 1))
  )
}

legends <- cowplot::plot_grid(
  get_leg(study),
  get_leg(disease),
  get_leg(dataset),
  get_leg(platform),
  get_leg(sex),
  get_leg(eth),
  get_leg(medication),
  get_leg(age + guides(fill = guide_colorbar(direction = "horizontal"))),
  get_leg(bmi + guides(fill = guide_colorbar(direction = "horizontal"))),
  get_leg(hba + guides(fill = guide_colorbar(direction = "horizontal"))),
  ncol = 1
)
# save --------------------------------------------------------------------

combined <- (study / disease / hba / bmi / age / sex / medication / eth / dataset / platform) &
  theme(
    plot.margin = margin(1, 1, 1, 1),
    panel.spacing = unit(0, "lines")
  )

ggplot2::ggsave(
  filename = here::here("islet_cartography_scrna/data/annotate/plot/meta_heat.pdf"),
  plot = combined,
  dpi = 300,
  width = 8,
  height = 3
)

ggplot2::ggsave(
  here::here("islet_cartography_scrna/data/annotate/plot/meta_heat_legends.pdf"),
  legends,
  dpi = 300,
  width = 12,
  height = 3
)


# barplots and histograms -------------------------------------------------
library(patchwork)

meta_donor <- meta_donor |> 
  dplyr::select(-ic_id_platform_adjusted_sample, -name, -ic_id_dataset, -ic_id_study, -library_prep, -cell_nuclei) |> 
  dplyr::distinct()

## ---- Age ----
p_age <- meta_donor |>  
  ggplot(aes(x = age_years)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 30,
                 fill = "grey",
                 color = "white") +
  geom_boxplot(aes(y = -2, group = disease_harmonized, fill = disease_harmonized),
               width = 3, outlier.size = 0.5) +
  scale_fill_manual(values = disease_color) +
  labs(x = "Age (years)", y = "Frequency", title = "Age") +
  my_theme() +
  theme(legend.position = "none")

## ---- BMI ----
p_bmi <- meta_donor |>  
  ggplot(aes(x = bmi)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 30,
                 fill = "grey",
                 color = "white") +
  geom_boxplot(aes(y = -2, group = disease_harmonized, fill = disease_harmonized),
               width = 3, outlier.size = 0.5) +
  scale_fill_manual(values = disease_color) +
  labs(x = "Body Mass Index (BMI)", y = "Frequency", title = "BMI") +
  my_theme() +
  theme(legend.position = "none")

## ---- HbA1c ----
p_hba <- meta_donor |>  
  ggplot(aes(x = hba_1_c_percent)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 30,
                 fill = "grey",
                 color = "white") +
  geom_boxplot(aes(y = -2, group = disease_harmonized, fill = disease_harmonized),
               width = 3, outlier.size = 0.5) +
  scale_fill_manual(values = disease_color) +
  labs(x = "HbA1c (%)", y = "Frequency", title = "HbA1c") +
  my_theme() +
  theme(legend.position = "none")

## ---- Disease ----
p_disease <- meta_donor |>
  count(disease_harmonized) |>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(disease_harmonized, perc, fill = disease_harmonized)) +
  geom_col() +
  scale_fill_manual(values = disease_color) +
  labs(x = "Disease", y = "Proportion of donors", title = "Disease") +
  my_theme() +
  theme(legend.position = "none")

## ---- Sex ----
p_sex <- meta_donor |>
  count(sex_predicted, disease_harmonized) |>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(disease_harmonized, perc, fill = sex_predicted)) +
  scale_fill_manual(values = c("#C2563A", "#3F7F93")) +
  geom_col(position = "dodge") +
  labs(x = "Disease", y = "Proportion of donors", title = "Sex") +
  my_theme() +
  theme(legend.position = "none")

## ---- Ethnicity ----
p_eth <- meta_donor |>
  count(ethnicity_broad_harmonized, disease_harmonized) |>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(disease_harmonized, perc, fill = ethnicity_broad_harmonized)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cols5) +
  labs(x = "Disease", y = "Proportionof donors", title = "Ethnicity") +
  my_theme() +
  theme(legend.position = "none")

  
## Combine (clean layout)

combined <- p_disease | p_sex | p_eth | p_age | p_bmi | p_hba

## save
ggplot2::ggsave(
  here::here("islet_cartography_scrna/data/annotate/plot/meta_histogram.pdf"),
  combined,
  width = 8,
  height = 2
)


# annotation --------------------------------------------------------------
# -----------------------------
# Disease composition
# -----------------------------
p_disease <- meta |>
  count(disease_harmonized, manual_annotation) |>
  ggplot(aes(x = disease_harmonized, y = n, fill = manual_annotation)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = manual_anno_colors) +
  labs(x = "Disease", y = "Proportion", fill = "Cell type") +
  my_theme()

# -----------------------------
# Dataset composition
# -----------------------------
p_dataset <- meta |>
  count(ic_id_dataset, manual_annotation) |>
  ggplot(aes(x = ic_id_dataset, y = n, fill = manual_annotation)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = manual_anno_colors) +
  labs(x = "Dataset", y = "Proportion", fill = "Cell type") +
  my_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

# -----------------------------
# Combine side-by-side
# -----------------------------
combined <- p_disease | p_dataset

# -----------------------------
# Save
# -----------------------------
ggplot2::ggsave(
  filename = here::here("islet_cartography_scrna/data/annotate/plot/celltype_composition.pdf"),
  plot = combined,
  dpi = 300,
  width = 8,
  height = 3
)

