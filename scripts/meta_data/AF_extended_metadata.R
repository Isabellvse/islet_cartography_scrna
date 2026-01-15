# Description -------------------------------------------------------------
# Here I will extract and harmonize medical information for donors where this is available
# as well as patch data
# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)


# Load --------------------------------------------------------------------
## The roughly cleaned metadata ----
# From scripts/meta_data/AB_combine_metadata.R
meta_list <- qs2::qs_read(here::here("islet_cartography_scrna/data/metadata/meta_list_1.qs2"))

## HPAP medical data ----
hpap_medications <- readxl::read_excel(here::here("islet_cartography_scrna/data_raw/meta/HPAP/Donor_Summary_176.xlsx"), 
                                           sheet = "medical_history_medications")  |> 
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(donor = donor_id)


# Preprocess Motakis ------------------------------------------------------
## Motakis ----
# Raw motakis data to match donor_ids (becuase in the GEO they are isletXXX, but in the supplementary data they are something different)
motakis <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Motakis/GSE221156_series_matrix.txt"), getGPL = F)
motakis <- Biobase::phenoData(motakis)@data

# Make motakis geo data fit with the supplementary data
motakis_test <- motakis |> 
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::select(sample = title,
                sex = sex_ch_1,
                age_years = age_ch_1,
                bmi_kg_m_2 = bmi_ch_1,
                "hb_a_1_c" = hba_1_c_ch_1,
                ancestry = race_ch_1,
                cause_of_death = cause_of_death_ch_1,
                unos = characteristics_ch_1_9, 
                run = characteristics_ch_1_12) |> 
  dplyr::mutate(unos = stringr::str_remove(unos, "unos: "),
                dplyr::across(c("age_years", "bmi_kg_m_2", "hb_a_1_c"), ~ base::as.numeric(.)), 
                ancestry = dplyr::case_when(ancestry == "White" ~ "E",
                                            ancestry == "AfricanAmerican" ~ "AA",
                                            ancestry == "Hispanic" ~ "H"))  |> 
  dplyr::filter(!grepl("/", sample)) |> 
  dplyr::select(sample, unos, "hb_a_1_c", "sex", "ancestry", "age_years")

# Motakis supplementary table 1
motakis_medication <- readxl::read_excel(here::here("islet_cartography_scrna/data_raw/meta/Motakis/tabel_1.xlsx"), skip = 1) |> 
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(donor = donor_id) |> dplyr::select(donor, unos, "hb_a_1_c", "sex", "ancestry", "age_years", diabetes_medication = diabetic_medication)

# combine the two tables
motakis_meds <- dplyr::right_join(x = motakis_medication, y = motakis_test) |> 
  dplyr::select(sample, diabetes_medication) |> 
  dplyr::mutate(diabetes_medication = dplyr::case_when(grepl("Metformin|oral", diabetes_medication) ~ "yes",
                                        .default = "no"))


# Preprocess HPAP medication ----------------------------------------------

## HPAP medication ----
hpap_donors <- c(meta_list[["HPAP_10x"]]$donor, 
                 meta_list[["HPAP_fluidigm"]]$donor, 
                 meta_list[["HPAP_patch_22"]]$donor, 
                 meta_list[["HPAP_patch_23"]]$donor) |> 
  base::unique()

# Medication the donors got for diabetes
# if medication was not know it was denoted as "no"
hpap_medications_harmonized <- hpap_medications |> 
  dplyr::rename(diabetes_medication = medication) |> 
  dplyr::filter(donor %in% hpap_donors & condition == "Type of Diabetes") |> 
  dplyr::mutate(diabetes_medication_harmonized = dplyr::case_when(
    BiocGenerics::grepl("Unknown", diabetes_medication) ~ "no",
    BiocGenerics::grepl("Insulin|Glipizide|Metformin|Sitaglipti|Sitaglipti|Lantus|mg unknown|", diabetes_medication) ~ "yes",
    .default = "no")) |> 
  dplyr::select(-compliance, - condition) 


# Add medication data to meta data ----------------------------------------
# Here I will only process those who contain information about medication
meta_list_extended <- meta_list |> 
  purrr::map_at(.at = "Lawlor", ~ dplyr::select(.x, name, donor, sample, donor_study, diabetes_medication = on_diabetes_medication)) |> 
  purrr::map_at(.at = "Son", ~ dplyr::select(.x, name, donor, sample, donor_study, diabetes_medication = treatment) |> 
                  dplyr::mutate(diabetes_medication_harmonized = dplyr::case_when(
                    BiocGenerics::grepl("Insulin|Metformin", diabetes_medication) ~ "yes",
                    .default = "no"))) |> 
  purrr::map_at(.at = "Wang_Sander", ~ dplyr::select(.x, name, donor, sample, donor_study, diabetes_medication = treatment_history) |> 
                  dplyr::mutate(diabetes_medication_harmonized = dplyr::case_when(
                    BiocGenerics::grepl("oral medication|meds", diabetes_medication) ~ "yes",
                    .default = "no"))) |> 
  purrr::map_at(.at = "Motakis", ~ dplyr::left_join(x = .x, y = motakis_meds) |> 
                  dplyr::select(name, donor, sample, donor_study, diabetes_medication)) |> 
  purrr::map_at(.at = "HPAP_10x", ~ dplyr::left_join(x = .x, y = hpap_medications_harmonized) |> 
                  dplyr::select(name, donor, sample, donor_study, diabetes_medication)) |> 
  purrr::map_at(.at = "HPAP_fluidigm", ~ dplyr::left_join(x = .x, y = hpap_medications_harmonized) |> 
                  dplyr::select(name, donor, sample, donor_study, diabetes_medication)) |> 
  purrr::map_at(.at = "HPAP_patch_22", ~ dplyr::left_join(x = .x, y = hpap_medications_harmonized) |> 
                  dplyr::select(name, donor, sample, donor_study, diabetes_medication)) |> 
  purrr::map_at(.at = "HPAP_patch_23", ~ dplyr::left_join(x = .x, y = hpap_medications_harmonized) |> 
                  dplyr::select(name, donor, sample, donor_study, diabetes_medication))


# Add patch-seq data to meta data -----------------------------------------
camunas <- c(
  "calcium_integral_normalizedto_cell_size_p_c_p_f",
  "calcium_integralat_first_depolarization_p_c",
  "capacitance_normalizedto_calcium_f_f_p_c",
  "cell_size_p_f",
  "cell_type_estimate_patching",
  "early_peak_calcium_current_amplitudeat_10_m_v_p_a",
  "first_depolarization_capacitance_f_f",
  "glucose_m_m",
  "half_inactivation_ca_current_m_v",
  "half_inactivation_sodium_current_m_v",
  "hyperpolarizationactivatedcurrent_at_140_m_v_p_a",
  "late_calcium_channel_conductance_p_s",
  "late_calcium_current_amplitudeat_10_m_v_p_a",
  "late_depolarization_capacitance",
  "normalized_early_peak_calcium_current_amplitude_p_a_p_f",
  "normalized_first_depolarization_capacitance_f_f_p_f",
  "normalized_hyperpolarizationactivatedcurrent_at_140_m_v_p_a_p_f",
  "normalized_late_calcium_channel_conductance_p_s_p_f",
  "normalized_late_calcium_current_amplitude_p_a_p_f",
  "normalized_late_depolarization_capacitance",
  "normalized_peak_sodium_current_amplitude_p_a_p_f",
  "normalized_sodium_channel_conductance_p_s_p_f",
  "normalized_total_capacitance_f_f_p_f",
  "peak_sodium_current_amplitudeat_0_m_v_p_a",
  "reversal_potentialby_cacurrents_m_v",
  "reversal_potentialbyramp_m_v",
  "sodium_channel_conductance_p_s",
  "solutionin_dish",
  "timefrom_dispersion_days",
  "timein_dish_min",
  "total_capacitance_f_f",
  "voltage_sodium_peak_current_m_v",
  "total_ieq",
  "ieq_per_pancreas_weight_ieq_g",
  "dna_content_ug",
  "insulin_content_ug",
  "insulin_content_per_ieq_ng_ieq",
  "insulin_to_dna_ratio"
)

# Vector 2 (to be corrected)
dai <- c(
  "glucose_m_m",
  "sodium_channel_conductance_p_s",
  "normalized_sodium_channel_conductance_p_s_p_f",
  "late_calcium_channel_conductance_p_s" = "late_ca_channel_conductance_p_s",
  "normalized_late_calcium_channel_conductance_p_s_p_f" = "normalized_late_ca_channel_conductance_p_s_p_f",
  "half_inactivation_ca_current_m_v",
  "total_capacitance_f_f",
  "normalized_total_capacitance_f_f_p_f",
  "first_depolarization_capacitance_f_f",
  "normalized_first_depolarization_capacitance_f_f_p_f",
  "late_depolarization_capacitance",
  "normalized_late_depolarization_capacitance",
  "calcium_integralat_first_depolarization_p_c",
  "calcium_integral_normalizedto_cell_size_p_c_p_f",
  "capacitance_normalizedto_calcium_f_f_p_c",
  "peak_sodium_current_amplitudeat_10_m_v_p_a",
  "normalized_peak_sodium_current_amplitudeat_10_m_v_p_a_p_f",
  "half_inactivation_sodium_current_m_v" = "half_inactivationof_sodium_current_m_v",
  "voltage_sodium_peak_current_m_v" = "voltagefor_sodium_peak_current_m_v",
  "early_peak_calcium_current_amplitudeat_10_m_v_p_a",
  "normalized_early_peak_calcium_current_amplitudeat_10_m_v_p_a_p_f",
  "late_calcium_current_amplitudeat_10_m_v_p_a",
  "normalized_late_calcium_current_amplitudeat_10_m_v_p_a_p_f",
  "lva_60_m_v_p_a",
  "normalized_lva_60_m_v_p_a_p_f",
  "hva_20_m_v_p_a",
  "normalized_hva_20_m_v_p_a_p_f",
  "reversal_potentialbyramp_m_v",
  "hyperpolarizationactivatedcurrent_at_140_m_v_p_a",
  "normalized_hyperpolarizationactivatedcurrent_at_140_m_v_p_a_p_f"
)


meta_list_extended <- meta_list_extended |> 
  purrr::map_at(.at = "Camunas", ~ dplyr::select(.x, name, donor, sample, donor_study, tidyselect::all_of(camunas)) |> 
                  dplyr::mutate(dplyr::across(-c(name, donor, sample, donor_study), ~ suppressWarnings(as.numeric(.))))) |> 
  purrr::map_at(.at = "Dai", ~ dplyr::select(.x, name, donor, sample, donor_study,tidyselect::all_of(dai)) |> 
                  dplyr::mutate(dplyr::across(-c(name, donor, sample, donor_study), ~ suppressWarnings(as.numeric(.)))))


# Save meta data ----------------------------------------------------------
## Combine this metadata
meds <- c("Lawlor", "Son", "Wang_Sander", "Motakis", "HPAP_10x", "HPAP_fluidigm", "HPAP_patch_22", "HPAP_patch_23")
patch <- c("Camunas", "Dai")

meta_list_extended[meds] |> 
  purrr::list_rbind() |> 
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata_harmonized/meta_harmonized_medication.csv"), delim = ",", col_names = TRUE)

meta_list_extended[patch] |> 
  purrr::list_rbind() |> 
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata_harmonized/meta_harmonized_patch.csv"), delim = ",", col_names = TRUE)

