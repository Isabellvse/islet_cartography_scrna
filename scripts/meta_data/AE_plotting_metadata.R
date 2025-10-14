# Description -------------------------------------------------------------
# Plotting meta data from donors that are left after QC
# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load --------------------------------------------------------------------
# STAR quality control
star_quality <- qs2::qs_read(here::here("islet_cartography_scrna/data/quality_control/star_quality_raw.qs2")) |> 
  # Add whether the method is droplet or plate based, ensure library prep is lower case
  purrr::modify_depth(1, ~ dplyr::mutate(., library_prep = base::tolower(library_prep), 
                                         platform = dplyr::case_when(library_prep %in% droplet_based ~ "droplet",
                                                                     library_prep %in% plate_based ~ "plate",
                                                                     library_prep %in% plate_based_bc ~ "plate_barcode")) %>% 
                        dplyr::relocate(platform, .after = "library_prep")) |>  
  purrr::list_rbind() |> dplyr::select(ic_id)

# Meta data - quality filtered data 
meta <- vroom::vroom(here::here("islet_cartography_scrna/data/AB_combined_metadata.csv"))

# raw meta data
meta_raw <- vroom::vroom(here::here("islet_cartography_scrna/data/metadata_harmonized/metadata_combined_harmonized.csv")) |> 
  dplyr::mutate(ic_id = ic_id_sample) |> 
  dplyr::select(starts_with("ic_id"), platform, name, study) |> 
  dplyr::distinct()


# How many samples and donor originally -----------------------------------
meta_star <- dplyr::left_join(star_quality, meta_raw, by = "ic_id")

# Generate the study overall annotation
meta_star <- meta_star |> dplyr::mutate(
  name_2 = case_when(
    grepl("^HPAP", name) ~ "HPAP",
    TRUE ~ as.character(name)))

# Generate the donor overall 
meta_star <- meta_star |> 
  dplyr::select(name_2, ic_id_donor, study)  |> 
  dplyr::distinct()  |> 
  dplyr::group_by(name_2) |> 
  dplyr::summarise(
    study_id = base::paste0(unique(study), collapse = ""),  # collapse into "896"
    .groups = "drop"
  ) |> 
  dplyr::right_join(y=meta_star) |> 
  dplyr::mutate(donor_number = str_extract(ic_id_donor, "(?<=_)\\d+$"),
                ic_id_donor_overall = paste("ic", study_id, donor_number, sep = "_")) |> 
  dplyr::relocate(ic_id_donor_overall, .after = ic_id_sample) |> 
  dplyr::select(-name_2, -donor_number, -study_id)

# number of donors (ic_19_11 has no donor id since is was excluded before ic_id generation) - because the authors
# stated that it failed sequencing
meta_star |> dplyr::select(ic_id_donor_overall) |> dplyr::distinct() |> dplyr::pull(ic_id_donor_overall) |> length()

# number of samples - here we count donors as samples for plate based methods
meta_star |> dplyr::mutate(sample = dplyr::case_when(platform == "droplet" ~ as.character(ic_id_sample),
                                                     .default = as.character(ic_id_donor))) |> 
  dplyr::select(sample) |> dplyr::distinct() |> dplyr::pull(sample) |> length()

# number of samples after QC
meta |> dplyr::mutate(sample = dplyr::case_when(platform %in% c("droplet", "plate_barcode") ~ as.character(ic_id_sample),
                                                     .default = as.character(ic_id_donor))) |> 
  dplyr::select(sample) |> dplyr::distinct() |> dplyr::pull(sample) |> length()

# number of donors after QC
meta |> dplyr::select(ic_id_donor_overall) |> dplyr::distinct() |> dplyr::pull(ic_id_donor_overall) |> length()
# Plots -------------------------------------------------------------------
test <- meta |>
  dplyr::select(ic_id_study_overall,
                doi,
                bioproject, 
                bmi,
                pmid,
                ic_id_donor_overall,
                age_years,
                cause_of_death_broad_harmonized,
                disease,
                ethnicity_broad_harmonized,
                ethnicity_sub_harmonized,
                sex_predicted,
                hba_1_c_percent,
                library_prep
                ) |> 
  dplyr::distinct() |> 
  dplyr::group_by(ic_id_donor_overall) |> 
  dplyr::mutate(library_prep = base::paste0(BiocGenerics::unique(library_prep), 
                                            collapse = ",")) |> 
  dplyr::ungroup() |> 
  dplyr::distinct()

donor_df <- meta |>
  dplyr::select(ic_id_donor_overall,
                age_years,
                cause_of_death_broad_harmonized,
                disease,
                bmi,
                ethnicity_broad_harmonized,
                ethnicity_sub_harmonized,
                sex_predicted,
                hba_1_c_percent
  ) |> 
  dplyr::distinct() |> 
  dplyr::mutate(disease_harmonized = dplyr::case_when(base::grepl("t2d", disease) ~ "t2d",
                                                   .default = base::as.character(disease)))

# Plotting donors ---------------------------------------------------------
pdf("meta_data_with_legends.pdf", width = 2, height = 2)

## ---- Disease not harmonized ----
p1 <- donor_df |> 
  dplyr::group_by(disease) |> 
  dplyr::tally() |> 
  dplyr::mutate(perc = round((n / sum(n))*100, 2)) |> 
  ggplot(aes(x = disease, y = perc, fill = disease)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Disease",
       y = "% of donors",
       title = "Disease distribution") +
  scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$disease)))) +
  my_theme() +
  ggplot2::theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

legend1 <- cowplot::get_legend(p1)
print(p1 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend1)

## ---- Disease harmonized ----
p2 <- donor_df |> 
  dplyr::group_by(disease_harmonized) |> 
  dplyr::tally() |> 
  dplyr::mutate(perc = round((n / sum(n))*100, 2)) |> 
  ggplot(aes(x = disease_harmonized, y = perc, fill = disease_harmonized)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Disease",
       y = "% of donors",
       title = "Disease distribution (harmonized)") +
  scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$disease_harmonized)))) +
  my_theme()

legend2 <- cowplot::get_legend(p2)
print(p2 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend2)

## ---- Sex ----
p3 <- donor_df |>
  dplyr::group_by(sex_predicted, disease_harmonized) |>
  dplyr::tally() |>
  dplyr::ungroup() |>
  dplyr::mutate(perc = round((n / sum(n))*100, 2)) |>
  ggplot(aes(x = disease_harmonized, y = perc, fill = sex_predicted)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Sex",
       y = "% of donors",
       title = "Sex distribution (Predicted)") +
  scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$sex_predicted)))) +
  my_theme()

legend3 <- cowplot::get_legend(p3)
print(p3 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend3)

## ---- Ethnicity broad ----
p4 <- donor_df |>
  dplyr::group_by(ethnicity_broad_harmonized, disease_harmonized) |>
  dplyr::tally() |> 
  dplyr::ungroup() |> 
  dplyr::mutate(perc = round((n / sum(n))*100, 2)) |> 
  ggplot(aes(x = disease_harmonized, y = perc, fill = ethnicity_broad_harmonized)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Sex",
       y = "% of donors",
       title = "Ethnicity (broad) distribution") +
  scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$ethnicity_broad_harmonized)))) + 
  my_theme()

legend4 <- cowplot::get_legend(p4)
print(p4 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend4)

## ---- Ethnicity detailed ----
p5 <- donor_df |>
  dplyr::group_by(ethnicity_sub_harmonized, disease_harmonized) |>
  dplyr::tally() |> 
  dplyr::ungroup() |> 
  dplyr::mutate(perc = round((n / sum(n))*100, 2)) |> 
  ggplot(aes(x = disease_harmonized, y = perc, fill = ethnicity_sub_harmonized)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Sex",
       y = "% of donors",
       title = "Ethnicity (sub-categories) distribution") +
  scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$ethnicity_sub_harmonized)))) + 
  my_theme()

legend5 <- cowplot::get_legend(p5)
print(p5 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend5)

## ---- Age ----
p6 <- donor_df |>  
  ggplot(aes(x = age_years)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 30,
                 fill = khroma::color("light")(1),
                 color = "white") +
  geom_boxplot(aes(y = -2, group = disease_harmonized, fill = disease_harmonized),
               width = 3,
               outlier.size = 0.5) +
  scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$disease_harmonized))),
                    name = "Disease") +
  labs(x = "Age (Years)",
       y = "% of donors",
       title = "Age Distribution") +
  my_theme()

legend6 <- cowplot::get_legend(p6)
print(p6 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend6)

## ---- BMI ----
p7 <- donor_df |>  
  ggplot(aes(x = bmi)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 30,
                 fill = khroma::color("light")(1),
                 color = "white") +
  geom_boxplot(aes(y = -2, group = disease_harmonized, fill = disease_harmonized),
               width = 3,
               outlier.size = 0.5) +
  scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$disease_harmonized))),
                    name = "Disease") +
  labs(x = "Body Mass Index (BMI)",
       y = "% of donors",
       title = "BMI Distribution") +
  my_theme()

legend7 <- cowplot::get_legend(p7)
print(p7 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend7)

## ---- HbA1c ----
p8 <- donor_df |>  
  ggplot(aes(x = hba_1_c_percent)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                 bins = 30,
                 fill = khroma::color("light")(1),
                 color = "white") +
  geom_boxplot(aes(y = -2, group = disease_harmonized, fill = disease_harmonized),
               width = 3,
               outlier.size = 0.5) +
  scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$disease_harmonized))),
                    name = "Disease") +
  labs(x = "Hemoglobin A1C (HbA1c) (%)",
       y = "% of donors",
       title = "HbA1c Distribution") +
  my_theme()

legend8 <- cowplot::get_legend(p8)
print(p8 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend8)

## ---- Cause of death broad ----
p9 <- donor_df |>
  dplyr::group_by(cause_of_death_broad_harmonized, disease_harmonized) |>
  dplyr::tally() |> 
  dplyr::ungroup() |> 
  dplyr::mutate(perc = round((n / sum(n))*100, 2)) |> 
  ggplot(aes(x = disease_harmonized, y = perc, fill = cause_of_death_broad_harmonized)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Sex",
       y = "% of donors",
       title = "Cause of death (broad) distribution") +
  scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$cause_of_death_broad_harmonized)))) + 
  my_theme()

legend9 <- cowplot::get_legend(p9)
print(p9 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend9)

## Library prep ----
library_df <- meta |> 
  dplyr::select(ic_id_study, library_prep) |> 
  dplyr::distinct() |> 
  dplyr::mutate(based = dplyr::case_when(library_prep %in% plate_based ~ "plate",
                                         library_prep %in% droplet_based ~ "droplet",
                                         library_prep %in% plate_based_bc ~ "plate"))

variable_order <- library_df |> 
  dplyr::group_by(library_prep, based) |> 
  dplyr::tally() |> 
  dplyr::arrange(based, desc(n)) |> 
  pull(library_prep)

study_order <- library_df |> 
  dplyr::mutate(based = factor(based, levels = c("droplet", "plate"))) |> 
  dplyr::arrange(based, library_prep) |> 
  dplyr::pull(ic_id_study) |> 
  base::unique()

p10 <- library_df |> 
  dplyr::mutate(ic_id_study = factor(ic_id_study, levels = study_order),
                library_prep= factor(library_prep, levels = variable_order)) |> 
  ggplot2::ggplot(aes(x = ic_id_study, y = library_prep, fill = based)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::labs(x = "Library preperation method",
       y = "Library preperation method distribution",
       title = "Datasets") +
  ggplot2::scale_fill_manual(values = khroma::color("light")(length(unique(donor_df$cause_of_death_broad_harmonized)))) + 
  my_theme() +
  ggplot2::theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

legend10 <- cowplot::get_legend(p10)
print(p10 + theme(legend.position = "none"))
grid::grid.newpage()
grid::grid.draw(legend10)

dev.off()


# pca regressions ---------------------------------------------------------
pcr <- vroom::vroom(here::here("islet_cartography_scrna/data/integrate/files/pcr_permutation_results.csv"))
pcr
# Overview of missing data ------------------------------------------------

#  Order variables by % missing
donor_df_deatiled <- meta |>
  dplyr::select(ic_id_donor_overall,
                age_years,
                cause_of_death_broad_harmonized,
                disease,
                ethnicity_broad_harmonized,
                ethnicity_sub_harmonized,
                sex_predicted,
                hba_1_c_percent,
                islet_culture_medium,
                islet_culture_medium_glucose_milimolar,
                dissociation_method,
                islet_isolation_enzyme,
                islet_fresh_frozen,
                islet_allocation_facility,
                treatment_patch,
                treatment_facs,
                cell_nuclei,
                strandedness,
                library_prep,
                sequencing_run,
                count_molecule,
                doi,
                year_public,
                tissue,
                islet_center,
                library_layout,
                count_quantification,
                geo_accession
  ) |> 
  dplyr::distinct() |> 
  dplyr::mutate(disease_harmonized = dplyr::case_when(base::grepl("t2d", disease) ~ "t2d",
                                                      .default = base::as.character(disease)))

variable_order <- donor_df_deatiled |> 
  summarise(across(-ic_id_donor_overall, ~ mean(is.na(.) | . == ""))) |> 
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") |> 
  arrange(desc(pct_missing)) |> 
  pull(variable)

#  Order donors by disease
donor_order <- donor_df_deatiled |> 
  dplyr::mutate(disease_harmonized = factor(disease_harmonized, levels = c("nd","pre","t2d"))) |> 
  dplyr::arrange(disease_harmonized) |> 
  dplyr::pull(ic_id_donor_overall) |> 
  base::unique()

# Prepare long-format data
donor_long <- donor_df_deatiled |> 
  dplyr::mutate(row_id = factor(ic_id_donor_overall, levels = donor_order)) |>  # x-axis ordered by disease
  dplyr::mutate(across(-c(ic_id_donor_overall, disease_harmonized, row_id), as.character)) |> 
  tidyr::pivot_longer(
    cols = -c(ic_id_donor_overall, disease_harmonized, row_id),
    names_to = "variable",
    values_to = "value"
  ) |> 
  dplyr::mutate(
    missing = is.na(value) | value == "",
    # For coloring: if missing, use "Missing"; else use disease
    fill_group = ifelse(missing, "Missing", disease_harmonized),
    variable = factor(variable, levels = variable_order)
  )

# Plot
missing_data <- ggplot2::ggplot(donor_long, aes(x = row_id, y = variable, fill = fill_group)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_manual(
    values = c(
      "nd" = "#8ABBEF", 
      "pre" = "#FF9B7C", 
      "t2d" = "#F5E48F",
      "Missing" = "white"
    ),
    name = "Disease / Missing"
  ) +
  ggplot2::labs(
    x = "Donors ordered by disease",
    y = "Variables (ordered by most missing)",
    title = "Missing Data"
  ) +
  my_theme() +
  ggplot2::theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Permutation test --------------------------------------------------------
pdf("covariate.pdf", width = 2, height = 2)
pcr_permutation_results |> 
  dplyr::mutate(covariate = factor(covariate, levels = pcr_permutation_results |> 
                                     dplyr::arrange(observed_r2) |> pull(covariate))) |> 
  ggplot2::ggplot(aes(x = observed_r2, y = covariate, fill = covariate))+
  geom_bar(stat="identity") +
  scale_fill_manual(values = khroma::color("light")(length(pcr_permutation_results$covariate)),
                    name = "Covariate") +
  my_theme() +
  theme(legend.position = "None")

dev.off()

