# Description -------------------------------------------------------------
# Plotting meta data from donors left from second around of QC (per cluster)
# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

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

# Load --------------------------------------------------------------------
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/paper_figures/files/obs.csv"))

# Preprocess --------------------------------------------------------------
meta_donor <- meta |> 
  dplyr::select(name, ic_id_donor_overall, 
                ic_id_platform_adjusted_sample, 
                ic_id_study, ic_id_dataset, library_prep, 
                disease_harmonized,
                sex_predicted, age_years, bmi, hba_1_c_percent, 
                ethnicity_broad_harmonized,
                cell_nuclei, platform) |> 
  dplyr::distinct() |> 
  dplyr::mutate(ethnicity_broad_harmonized = dplyr::case_when(ethnicity_broad_harmonized == "of_european_descent" ~ "EU",
                                                              ethnicity_broad_harmonized == "of_african_descent" ~ "AF",
                                                              ethnicity_broad_harmonized == "of_latin_american_hispanic_descent" ~ "LAH",
                                                              ethnicity_broad_harmonized == "of_asian_descent" ~ "AS",
                                                              .default = NA),
                platform = dplyr::case_when(platform == "plate_barcode" ~ "plate",
                                            .default = as.character(platform)))


# Heatmap overview --------------------------------------------------------
## Donor order ----
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

## Generate barplots for each variable ----
study   <- tile_row(meta_donor, ic_id_study,   "Study") + scale_fill_manual(values = cols22)
disease <- tile_row(meta_donor, disease_harmonized, "Disease") + ggplot2::scale_fill_manual(values = disease_color)
sex <- tile_row(meta_donor, sex_predicted, "Sex") + scale_fill_manual(values = c("#C2563A", "#3F7F93"))
eth <- tile_row(meta_donor, ethnicity_broad_harmonized, "Ethnicity") + scale_fill_manual(values = cols5)
age  <- tile_row(meta_donor, age_years, "Age (years)") + ggplot2::scale_fill_gradient(low = "white", high = "#5D3FD3")
bmi  <- tile_row(meta_donor, bmi, "BMI") + ggplot2::scale_fill_gradient(low = "white", high = "#FF5F1F")
hba <- tile_row(meta_donor, hba_1_c_percent, "HbA1c (%)") + scale_fill_gradientn(
  colours = c("forestgreen", "white", "gold", "#ae0000"),
  values = scales::rescale(c(4.5, 5.7, 6.5, 13.1)),
  breaks = c(4.5, 5.7, 6.5, 13.1),
  name = NULL)

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

cell_nuc <- meta_donor |>
  count(ic_id_donor_overall, cell_nuclei) |>
  mutate(ic_id_donor_overall = factor(ic_id_donor_overall, levels = donor_order)) |>
  ggplot(aes(x = ic_id_donor_overall, y = n, fill = cell_nuclei)) +
  geom_col(position = "fill", width = 1, color = NA, linewidth = 0.3) +
  labs(y = "Material", fill = NULL) +
  scale_fill_manual(values = cols5) +
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
  count(ic_id_donor_overall, platform) |>
  mutate(ic_id_donor_overall = factor(ic_id_donor_overall, levels = donor_order)) |>
  ggplot(aes(x = ic_id_donor_overall, y = n, fill = platform)) +
  geom_col(position = "fill", width = 1, color = NA, linewidth = 0.3) +
  labs(y = "Platform", fill = NULL) +
  scale_fill_manual(values = cols5) +
  my_theme() +
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  )


## Combine barplots --------------------------------------------------------
combined <- (study / dataset / disease / hba / bmi / age / sex  / eth / platform / cell_nuc) &
  theme(
    plot.margin = margin(1, 1, 1, 1),
    panel.spacing = unit(0, "lines")
  )


# Legends -----------------------------------------------------------------
get_leg <- function(p) {
  cowplot::get_legend(
    p +
      theme(legend.position = "left",
            legend.text = element_text(angle = 90, hjust = 0.5, vjust = 0.5)))
}

legends <- cowplot::plot_grid(
  get_leg(disease),
  get_leg(hba),
  get_leg(bmi),
  get_leg(age),
  get_leg(sex),
  get_leg(eth),
  get_leg(platform),
  get_leg(cell_nuc),
  nrow = 1
)


# Save plots --------------------------------------------------------------
ggplot2::ggsave(
  filename = here::here("islet_cartography_scrna/data/paper_figures/plot/overview_heatmap.pdf"),
  plot = combined,
  dpi = 300,
  width = 100,
  units = "mm",
  height = 40
)

ggplot2::ggsave(
  here::here("islet_cartography_scrna/data/paper_figures/plot/overview_heatmap_legends.pdf"),
  legends,
  dpi = 300,
  width = 100,
  units = "mm",
  height = 40
)



# Barplots ----------------------------------------------------------------
# only keep donor specific data now
meta_donor <- meta_donor |> 
  dplyr::select(-ic_id_platform_adjusted_sample, -name, -ic_id_dataset, -ic_id_study, -library_prep, -cell_nuclei) |> 
  dplyr::distinct()

## ---- Age ----
p_age <- meta_donor |>  
  ggplot(aes(x = age_years)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 30,
                 fill = "grey") +
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
                 fill = "grey") +
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
                 fill = "grey") +
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

combined <- (p_disease | p_sex | p_eth) / (p_age | p_bmi | p_hba)

## save
ggplot2::ggsave(
  here::here("islet_cartography_scrna/data/paper_figures/plot/meta_histogram.pdf"),
  combined,
  units = "mm",
  width = 80,
  height = 80
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


