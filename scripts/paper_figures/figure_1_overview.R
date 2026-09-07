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
iridescent <- khroma::color("iridescent")
# Load --------------------------------------------------------------------
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/obs.csv")) 

# Statistical analysis ----------------------------------------------------
df <- meta |> 
  dplyr::select(ic_id_donor_overall, bmi, hba_1_c_percent, disease_hba1c) |> 
  dplyr::distinct() |> 
  dplyr::mutate(disease_hba1c = factor(disease_hba1c, levels = c("nd", "pre", "t2d")))

# summarise data
group_by(df, disease_hba1c) %>%
  summarise(
    count = n(),
    mean = mean(hba_1_c_percent, na.rm = TRUE),
    sd = sd(hba_1_c_percent, na.rm = TRUE),
    median = median(hba_1_c_percent, na.rm = TRUE),
    IQR = IQR(hba_1_c_percent, na.rm = TRUE)
  )
hba_data <- df |> 
  dplyr::select(-bmi) |> 
  tidyr::drop_na()

df |> 
  tidyr::pivot_longer(tidyselect::where(is.numeric)) |> 
  tidyr::nest(.by = "name") |> 
  dplyr::mutate(overall = purrr::map(data, \(df){
    df |> 
      rstatix::kruskal_test(value ~ disease_hba1c) |> 
      dplyr::left_join(df |> 
                         rstatix::kruskal_effsize(value ~ disease_hba1c) |> 
                         dplyr::select(-method)) |> 
      dplyr::select(-".y.")
  }),
  pairwise = purrr::map(data, \(df){
    df |> 
      rstatix::wilcox_test(value ~ disease_hba1c) |> 
      dplyr::mutate(method = "wilcox") |> 
      dplyr::left_join(df |> 
                         rstatix::wilcox_effsize(value ~ disease_hba1c)) |> 
      dplyr::select(-".y.")
  })) |> 
  dplyr::select(-data) |> 
  tidyr::pivot_longer(c(overall, pairwise), names_to = "test_type", values_to = "result") |> 
  tidyr::unnest(result) |> 
  vroom::vroom_write(here::here("islet_cartography_scrna/data/blood_and_bmi_statistical_test.csv"))


# Preprocess --------------------------------------------------------------
meta_donor <- meta |> 
  dplyr::select(name, ic_id_donor_overall, 
                ic_id_platform_adjusted_sample, 
                ic_id_study, ic_id_dataset, library_prep, 
                disease_hba1c,
                sex_predicted, age_years, bmi, hba_1_c_percent, 
                ethnicity_broad_harmonized,
                cell_nuclei, platform) |> 
  dplyr::distinct() |> 
  dplyr::mutate(ethnicity_broad_harmonized = dplyr::case_when(ethnicity_broad_harmonized == "of_european_descent" ~ "EU",
                                                              ethnicity_broad_harmonized == "of_african_descent" ~ "AF",
                                                              ethnicity_broad_harmonized == "of_latin_american_hispanic_descent" ~ "LH",
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
    ic_id_dataset,
    disease_hba1c,
    hba_1_c_percent,
    bmi,
    age_years,
    sex_predicted
  ) |>
  dplyr::distinct(ic_id_donor_overall, .keep_all = TRUE) |>
  dplyr::pull(ic_id_donor_overall)

## Generate barplots for each variable ----
study   <- tile_row(meta_donor, ic_id_study,   "Study") + scale_fill_manual(values = cols22)
disease <- tile_row(meta_donor, disease_hba1c, "Disease") + ggplot2::scale_fill_manual(values = disease_color)
sex <- tile_row(meta_donor, sex_predicted, "Sex") + scale_fill_manual(values = c("#C2563A", "#3F7F93"))
eth <- tile_row(meta_donor, ethnicity_broad_harmonized, "Ethnicity") + scale_fill_manual(values = cols5)
age  <- tile_row(meta_donor, age_years, "Age (years)") + ggplot2::scale_fill_gradient(low = "white", high = "#5D3FD3")
bmi  <- tile_row(meta_donor, bmi, "BMI") + ggplot2::scale_fill_gradient(low = "white", high = "#FF5F1F")
hba <- tile_row(meta_donor, hba_1_c_percent, "HbA1c (%)") + ggplot2::scale_fill_gradient(low = "white", high = "#FF5F1F")

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
  geom_boxplot(aes(y = -2, group = disease_hba1c, fill = disease_hba1c),
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
  geom_boxplot(aes(y = -2, group = disease_hba1c, fill = disease_hba1c),
               width = 3, outlier.size = 0.5) +
  scale_fill_manual(values = disease_color) +
  labs(x = "Body Mass Index (BMI)", y = "Frequency", title = "BMI") +
  my_theme() +
  theme(legend.position = "none", 
        axis.title.y = ggplot2::element_blank())

## ---- HbA1c ----
p_hba <- meta_donor |>  
  ggplot(aes(x = hba_1_c_percent)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 30,
                 fill = "grey") +
  geom_boxplot(aes(y = -2, group = disease_hba1c, fill = disease_hba1c),
               width = 3, outlier.size = 0.5) +
  scale_fill_manual(values = disease_color) +
  labs(x = "HbA1c (%)", y = "Frequency", title = "HbA1c") +
  my_theme() +
  theme(legend.position = "none",
        axis.title.y = ggplot2::element_blank())

## ---- Disease ----
p_disease <- meta_donor |>
  count(disease_hba1c) |>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(disease_hba1c, perc, fill = disease_hba1c)) +
  geom_col() +
  scale_fill_manual(values = disease_color) +
  labs(x = "Disease", y = "Proportion of donors", title = "Disease") +
  my_theme() +
  theme(legend.position = "none")

## ---- Sex ----
p_sex <- meta_donor |>
  count(sex_predicted, disease_hba1c) |>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(disease_hba1c, perc, fill = sex_predicted)) +
  scale_fill_manual(values = c("#C2563A", "#3F7F93")) +
  geom_col(position = "dodge") +
  labs(x = "Disease", y = "Proportion of donors", title = "Sex") +
  my_theme() +
  theme(legend.position = "none",
        axis.title.y = ggplot2::element_blank())

## ---- Ethnicity ----
p_eth <- meta_donor |>
  count(ethnicity_broad_harmonized, disease_hba1c) |>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(disease_hba1c, perc, fill = ethnicity_broad_harmonized)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cols5) +
  labs(x = "Disease", y = "Proportionof donors", title = "Ethnicity") +
  my_theme() +
  theme(legend.position = "none",
        axis.title.y = ggplot2::element_blank())


## Combine (clean layout)

combined <- (p_disease | p_sex | p_eth) / (p_age | p_bmi | p_hba) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1), guides = "collect") &
  theme(plot.margin = margin(0, 0, 0, 0))
combined

## save
ggplot2::ggsave(
  here::here("islet_cartography_scrna/data/paper_figures/plot/meta_histogram.pdf"),
  combined,
  units = "mm",
  width = 80,
  height = 70
)


# broad cell type per donor -----------------------------------------------------
p1 <- meta |> 
  dplyr::group_by(cell_type_broad) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(total = sum(n),
                pct = n / total) |> 
  ggplot2::ggplot(ggplot2::aes(x = "x", y = pct, fill = cell_type_broad)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = cell_type_broad_colors) +
  my_theme() +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank())

n_cells_donor <- meta |> 
  dplyr::group_by(ic_id_donor_overall) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::ungroup() |> 
  dplyr::arrange(dplyr::desc(n)) |> 
  dplyr::mutate(rank = dplyr::row_number())

p2 <- n_cells_donor |> 
  ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(ic_id_donor_overall, rank), y = n)) +
  ggplot2::geom_bar(stat = "identity", width = 1) +
  ggplot2::labs(y = "Total cells") +
  my_theme() +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank())

p3 <- meta |> 
  dplyr::group_by(cell_type_broad, ic_id_donor_overall) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(ic_id_donor_overall) |> 
  dplyr::mutate(total = sum(n),
                pct = n / total) |> 
  dplyr::ungroup() |> 
  ggplot2::ggplot(ggplot2::aes(x = forcats::fct_relevel(ic_id_donor_overall, 
                                                        n_cells_donor$ic_id_donor_overall), 
                               y = pct, fill = cell_type_broad)) +
  ggplot2::geom_bar(stat = "identity", width = 1) +
  ggplot2::scale_fill_manual(values = cell_type_broad_colors) +
  my_theme() +
  ggplot2::theme( 
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank()) 

pdf(here::here("islet_cartography_scrna/data/paper_figures/plot/celltype_broad_composition_per_donor.pdf"),
    width = 5,
    height = 2)
((patchwork::plot_spacer() / p1) + patchwork::plot_layout(heights = c(0.3, 1)) |
    (p2 / p3) + patchwork::plot_layout(heights = c(0.3, 1))) +
  patchwork::plot_layout(widths = c(0.1, 1), guides = "collect")
dev.off()

# check that axis order is the same
identical(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$x$get_labels(),
ggplot2::ggplot_build(p3)$layout$panel_params[[1]]$x$get_labels())


# cell type annotation per donor ------------------------------------------
## Overall percentage ----
p1 <- meta |> 
  dplyr::group_by(cell_type) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(total = sum(n),
                pct = n / total) |> 
  ggplot2::ggplot(ggplot2::aes(x = "x", y = pct, fill = cell_type)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = cell_type_colors) +
  my_theme() +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank())

n_cells_donor <- meta |> 
  dplyr::group_by(ic_id_donor_overall) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::ungroup() |> 
  dplyr::arrange(dplyr::desc(n)) |> 
  dplyr::mutate(rank = dplyr::row_number())

p2 <- n_cells_donor |> 
  ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(ic_id_donor_overall, rank), y = n)) +
  ggplot2::geom_bar(stat = "identity", width = 1) +
  ggplot2::labs(y = "Total cells") +
  my_theme() +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank())

p3 <- meta |> 
  dplyr::group_by(cell_type, ic_id_donor_overall) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(ic_id_donor_overall) |> 
  dplyr::mutate(total = sum(n),
                pct = n / total) |> 
  dplyr::ungroup() |> 
  ggplot2::ggplot(ggplot2::aes(x = forcats::fct_relevel(ic_id_donor_overall, 
                                                        n_cells_donor$ic_id_donor_overall), 
                               y = pct, fill = cell_type)) +
  ggplot2::geom_bar(stat = "identity", width = 1) +
  ggplot2::scale_fill_manual(values = cell_type_colors) +
  my_theme() +
  ggplot2::theme( 
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()) 

pdf(here::here("islet_cartography_scrna/data/paper_figures/plot/celltype_composition_per_donor.pdf"),
    width = 5,
    height = 2)
((patchwork::plot_spacer() / p1) + patchwork::plot_layout(heights = c(0.3, 1)) |
    (p2 / p3) + patchwork::plot_layout(heights = c(0.3, 1))) +
  patchwork::plot_layout(widths = c(0.1, 1), guides = "collect")
dev.off()

# check that axis order is the same
identical(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$x$get_labels(),
          ggplot2::ggplot_build(p3)$layout$panel_params[[1]]$x$get_labels())


# Number of cell per donor - dataset ------------------------------------------------
# rank dataset
dataset_order <- meta |> 
  dplyr::group_by(ic_id_dataset, ic_id_donor_overall) |> 
  dplyr::tally() |> 
  dplyr::group_by(ic_id_dataset) |>
  dplyr::summarise(max = median(n)) |> 
  dplyr::ungroup() |> 
  dplyr::arrange(max) |> 
  dplyr::pull(ic_id_dataset)

p1 <- meta |> 
  dplyr::group_by(ic_id_dataset, ic_id_donor_overall) |> 
  dplyr::tally() |> 
  dplyr::group_by(ic_id_dataset) |> 
  dplyr::mutate(n_donor = dplyr::n()) |> 
  dplyr::ungroup() |> 
  ggplot2::ggplot(ggplot2::aes(x = n, y = forcats::fct_relevel(ic_id_dataset, dataset_order))) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = n_donor), outlier.size = 0.5, 
                        median.linewidth = 0.5) +
  ggplot2::scale_x_log10(
    breaks = scales::trans_breaks("log10", \(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ggplot2::scale_fill_gradientn(colors = iridescent(5)) +
  ggplot2::annotation_logticks(sides = "bottom", outside = T, 
                               short = unit(0.05, "cm"), 
                               mid = unit(0.1, "cm"), 
                               long = unit(0.15, "cm"))  +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_upper(x) |> 
                              stringr::str_replace("_", " ")) +
  ggplot2::labs(y = "Dataset",
                x = "Cell per donor") +
  my_theme() +
  ggplot2::coord_cartesian(clip = "off")

p2 <- meta |> 
  dplyr::distinct(ic_id_dataset, ic_id_donor_overall) |> 
  dplyr::group_by(ic_id_dataset) |> 
  dplyr::tally(name = "n_donor") |> 
  dplyr::ungroup() |> 
  ggplot2::ggplot(ggplot2::aes(x = n_donor, y = forcats::fct_relevel(ic_id_dataset, dataset_order),
                               fill = n_donor)) + 
  ggplot2::geom_bar(stat="identity") +
  ggplot2::scale_fill_gradientn(colors = iridescent(5)) +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_upper(x) |> 
                              stringr::str_replace("_", " ")) +
  ggplot2::labs(y = "Dataset",
                x = "Cell per donor") +
  my_theme() +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())

pdf(here::here("islet_cartography_scrna/data/paper_figures/plot/cell_per_donor_dataset.pdf"),
    width = 3,
    height = 2)
p1 + p2 + patchwork::plot_layout(guides = "collect", widths = c(1, 0.3))
dev.off()


# Number of cell per donor study ------------------------------------------------
# rank dataset
dataset_order <- meta |> 
  dplyr::group_by(ic_id_study, ic_id_donor_overall) |> 
  dplyr::tally() |> 
  dplyr::group_by(ic_id_study) |>
  dplyr::summarise(max = median(n)) |> 
  dplyr::ungroup() |> 
  dplyr::arrange(max) |> 
  dplyr::pull(ic_id_study)

p1 <- meta |> 
  dplyr::group_by(ic_id_study, ic_id_donor_overall) |> 
  dplyr::tally() |> 
  dplyr::group_by(ic_id_study) |> 
  dplyr::mutate(n_donor = dplyr::n()) |> 
  dplyr::ungroup() |> 
  ggplot2::ggplot(ggplot2::aes(x = n, y = forcats::fct_relevel(ic_id_study, dataset_order))) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = n_donor), outlier.size = 0.5, 
                        median.linewidth = 0.5) +
  ggplot2::scale_x_log10(
    breaks = scales::trans_breaks("log10", \(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ggplot2::scale_fill_gradientn(colors = iridescent(5)) +
  ggplot2::annotation_logticks(sides = "bottom", outside = T, 
                               short = unit(0.05, "cm"), 
                               mid = unit(0.1, "cm"), 
                               long = unit(0.15, "cm"))  +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_upper(x) |> 
                              stringr::str_replace("_", " ")) +
  ggplot2::labs(y = "Study",
                x = "Cell per donor") +
  my_theme() +
  ggplot2::coord_cartesian(clip = "off")

p2 <- meta |> 
  dplyr::distinct(ic_id_study, ic_id_donor_overall) |> 
  dplyr::group_by(ic_id_study) |> 
  dplyr::tally(name = "n_donor") |> 
  dplyr::ungroup() |> 
  ggplot2::ggplot(ggplot2::aes(x = n_donor, y = forcats::fct_relevel(ic_id_study, dataset_order),
                               fill = n_donor)) + 
  ggplot2::geom_bar(stat="identity") +
  ggplot2::scale_fill_gradientn(colors = iridescent(5)) +
  ggplot2::scale_y_discrete(labels = function(x) stringr::str_to_upper(x) |> 
                              stringr::str_replace("_", " ")) +
  ggplot2::labs(y = "Dataset",
                x = "Cell per donor") +
  my_theme() +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())

pdf(here::here("islet_cartography_scrna/data/paper_figures/plot/cell_per_donor_study.pdf"),
    width = 3,
    height = 2)
p1 + p2 + patchwork::plot_layout(guides = "collect", widths = c(1, 0.3))
dev.off()

