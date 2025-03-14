# Description -------------------------------------------------------------
# Here I process all metadata files to the same format

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/metadata/"))

# Avrahami ----------------------------------------------------------------
# Load and process GEO data from two different files, renaming columns for clarity
a_1 <- GEOquery::getGEO(filename = here::here("islet_cartography_scrna/data_raw/meta/Avrahami/GSE154126-GPL11154_series_matrix.txt"), getGPL = F)
a_1 <- Biobase::phenoData(a_1)@data |>
  dplyr::rename(sample = title,
                donor = source_name_ch1,
                inferred_cell_type = `celltype:ch1`)

a_2 <- GEOquery::getGEO(filename = here::here("islet_cartography_scrna/data_raw/meta/Avrahami/GSE154126-GPL16791_series_matrix.txt"), getGPL = F)
a_2 <- Biobase::phenoData(a_2)@data |>
  dplyr::rename(sample = title,
                donor = source_name_ch1,
                inferred_cell_type = `celltype:ch1`)

# Combine the processed data and clean up the sample and cell type information
a_bind <- dplyr::bind_rows(a_1, a_2) |>
  dplyr::mutate(sample = stringr::str_replace(sample, ": ", "_"),
                inferred_cell_type = stringr::str_remove(inferred_cell_type, "pancreatic islets;"))

# Merge the supplementary data with the combined GEO data
a_meta <- readxl::read_xlsx(here::here("islet_cartography_scrna/data_raw/meta/Avrahami/supplementary_table_1.xlsx")) |>
  dplyr::filter(!grepl("\\*|Toddler|Newborn|Adolescent", Condition)) |>
  dplyr::rename_all(snakecase::to_snake_case) |>
  dplyr::rename(age_years = age,
                disease = condition,
                donor = donor_id) |>
  dplyr::left_join(y = a_bind) |>
  # Convert character "NA" values to true NA values and add new columns for further analysis
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::na_if(., "NA"))) |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::mutate(name = "Avrahami",
                age_years = base::as.numeric(stringr::str_remove(age_years, "y")),
                donor_study = base::paste0(name, "_", donor),
                bmi = base::as.numeric(bmi),
                ethnicity = snakecase::to_snake_case(ethnicity),
                gender = dplyr::case_when(gender=="male" ~ "m",
                                       gender == "female" ~ "f"),
                disease = dplyr::case_when(disease == "Adult" ~ "nd",
                                           disease == "T2D" ~ "t2d"),
                islet_center = "Integrated Islet Distribution Program (IIDP)|the Human Pancreas Analysis Program (HPAP) consortium") |>
  # Filter out specific donor samples and save the final data frame to a CSV file
  dplyr::filter(!donor == "* samples reported in our previous study (Wang, 2016)") |>
  dplyr::rename(species = organism_ch_1,
                molecule = molecule_ch_1,
                tissue = tissue_ch_1) |>
 dplyr::select(-dplyr::contains("data_processing"),
               -dplyr::contains("characteristics_"),
               -dplyr::contains("supplementary"),
               -taxid_ch_1,
               -data_row_count,
               -dplyr::contains("category"))

# save metadata
a_meta |> vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Avrahami.csv"), delim = ",", col_names = TRUE)

# save column classes
purrr::map(a_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Avrahami_class.csv"), delim = ",", col_names = TRUE)
# Baron -------------------------------------------------------------------
# Load and process GEO data from the specified file, extracting phenotype data
baron <- GEOquery::getGEO(filename = here::here("islet_cartography_scrna/data_raw/meta/Baron/GSE84133-GPL16791_series_matrix.txt"), getGPL = F)
baron <- Biobase::phenoData(baron)@data

# Convert relevant columns to numeric and standardize donor and sample names
baron_meta <- baron |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(sample = title,
         sex = "sex_ch_1",
         age_years = "age_ch_1",
         bmi = "bmi_ch_1",
         disease = type_2_diabetes_mellitus_ch_1) |>
  dplyr::mutate(sample = stringr::str_remove(string = sample, pattern = "human pancreatic islets, ") |> snakecase::to_snake_case(),
                donor = sample,
                age_years = base::as.numeric(age_years),
                bmi = base::as.numeric(bmi),
         # Adjust the sex column to a simplified format and add new columns for further analysis
         sex = dplyr::case_when(sex=="male" ~ "m",
                         sex == "female" ~ "f"),
         name = "Baron",
         donor_study =base::paste0(name, "_", donor),
         disease = dplyr::case_when(disease == "No" ~ "nd",
                                    disease == "Yes" ~ "t2d")) |>
  # Filter out samples with age less than 20 years and save the final data frame to a CSV file
  dplyr::filter(!age_years < 20) |>
  dplyr::rename(species = organism_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -taxid_ch_1,
                -data_row_count,
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Baron.csv"), delim = ",", col_names = TRUE)

purrr::map(baron_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Baron_class.csv"), delim = ",", col_names = TRUE)

# Camunas -----------------------------------------------------------------
# Load and process metadata from two different files, removing completely empty columns
camunas_facs <- utils::read.delim(here::here("islet_cartography_scrna/data_raw/meta/Camunas/patchclamp_FACS_human.metadata.tab")) |>
  # Remove completely empty columns
  dplyr::select(dplyr::where(~ !base::all(base::is.na(.)))) |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(sample = name,
         donor = donor_id,
         age_years = age,
         disease = diabetes_status,
         bmi = body_mass_index_bmi,
         "hba1c_percent" = glycated_hemoglobin_hb_a_1_c,
         serum_glucose = serum_glucose_measurements) |>
  # Convert character "NA" and empty strings to true NA values, and convert specified columns to numeric
  dplyr::mutate(dplyr::across(base::is.character, ~ dplyr::na_if(., "NA")),
                dplyr::across(base::is.character, ~ dplyr::na_if(., "")),
                across(c("age_years", "bmi", "hba1c_percent", "dna_content_ug"), ~ as.numeric(.)))

camunas_cry <- utils::read.delim(here::here("islet_cartography_scrna/data_raw/meta/Camunas/patchclamp_wcryo_human.metadata.tab")) |>
  # Remove completely empty columns
  dplyr::select(dplyr::where(~ !base::all(base::is.na(.)))) |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(sample = name,
                donor = donor_id,
                age_years = age,
                disease = diabetes_status,
                bmi = body_mass_index_bmi,
                "hba1c_percent" = glycated_hemoglobin_hb_a_1_c,
                serum_glucose = serum_glucose_measurements) |>
  # Convert character "NA" and empty strings to true NA values, and convert specified columns to numeric
  dplyr::mutate(dplyr::across(base::is.character, ~ dplyr::na_if(., "NA")),
                dplyr::across(base::is.character, ~ dplyr::na_if(., " ")),
                dplyr::across(c("age_years", "bmi", "hba1c_percent", "dna_content_ug"), ~ base::as.numeric(.)))

# Combine the processed data from both files and add new columns
camunas_meta <- camunas_facs |>
  dplyr::full_join(camunas_cry) |>
  dplyr::mutate(name = "Camunas",
                donor_study = paste0(name, "_", donor)) |>
  dplyr::relocate("name", donor) |>
  # Filter out specific conditions and age groups, and standardize disease and sex values
  dplyr::filter(!disease == "cryo_T1D",
                !age_years < 20) |>
  # Renaming for standardization
  dplyr::mutate(disease = dplyr::case_when(disease %in% c("cryo_healthy", "healthy") ~ "nd",
                                           disease == "T2D" ~ "t2d",
                                           disease == "prediabetic" ~ "pre",
                                           disease == "T2Dreversed" ~ "t2d_reverse",
                                           disease == "elevated HbA1c" ~ "elevated_hba1c"),
                sex = snakecase::to_snake_case(sex)) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Camunas.csv"), delim = ",", col_names = TRUE)

purrr::map(camunas_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Camunas_class.csv"), delim = ",", col_names = TRUE)

# Enge --------------------------------------------------------------------
enge <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Enge/GSE81547_series_matrix.txt"), getGPL = F)
enge <- Biobase::phenoData(enge)@data

enge_meta <- enge |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(sample = title,
         gender = "gender_ch_1",
         age_years = "donor_age_ch_1",
         inferred_cell_type = "inferred_cell_type_ch_1") |>
  # Adding information from Supplementary Data Table 1
  dplyr::mutate(gender = dplyr::case_when(gender == "male" ~ "m",
                            gender == "female" ~ "f"),
         donor = stringr::str_remove(sample, "_[^_]*$"),
         name = "Enge",
         disease = "nd",
         bmi = dplyr::case_when(donor == "21yr_male" ~ 28.4,
                         donor == "38yr_female" ~ 29.5,
                         donor == "44yr_female" ~ 23.8,
                         donor == "54yr_male" ~ 27.9,
                         donor == "22yr_male" ~ 24.8),
         ethnicity = dplyr::case_when(donor == "21yr_male" ~ "caucasian",
                               donor == "38yr_female" ~ "caucasian",
                               donor == "44yr_female" ~ "african_american",
                               donor == "54yr_male" ~ "caucasian",
                               donor == "22yr_male" ~ "asian"),
         cause_of_death = dplyr::case_when(donor == "21yr_male" ~ "head_trauma",
                                    donor == "38yr_female" ~ "stroke",
                                    donor == "44yr_female" ~ "stroke",
                                    donor == "54yr_male" ~ "anorexia",
                                    donor == "22yr_male" ~ "head_trauma"),
         donor_study = base::paste0(name, "_", donor),
         age_years = base::as.numeric(age_years),
         islet_center = "Integrated Islet Distribution Network (IIDP)|National Diabetes Research Institute (NDRI)|UCSF Islet Isolation Core (San Francisco, CA USA)|International Institute for the Advancement of Medicine
(IIAM)") |>
  dplyr::filter(age_years > 20) |>
  dplyr::rename(tissue = source_name_ch_1,
                species = organism_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -taxid_ch_1,
                -data_row_count,
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Enge.csv"), delim = ",", col_names = TRUE)

purrr::map(enge_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Enge_class.csv"), delim = ",", col_names = TRUE)

# Fang --------------------------------------------------------------------
fang <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Fang/GSE101207-GPL11154_series_matrix.txt"), getGPL = F)
fang <- Biobase::phenoData(fang)@data

fang_meta <- fang |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(sample = title,
         gender = gender_ch_1,
         age_years = donor_age_ch_1,
         bmi = bmi_ch_1,
         disease = source_name_ch_1) |>
  dplyr::mutate(name = "Fang",
         donor = sample,
         donor_study = base::paste0(name, "_", donor),
         gender = snakecase::to_snake_case(gender),
         disease = dplyr::case_when(base::grepl("Non-diabetic", disease) ~ "nd",
                                    base::grepl("Diabetic", disease) ~ "t2d"),
         dplyr::across(c("bmi", "age_years"), ~ base::as.numeric(.)),
         islet_center = "Prodo Laboratories, Inc") |>
  dplyr::rename(tissue = tissue_ch_1,
                species = organism_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -taxid_ch_1,
                -data_row_count,
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Fang.csv"), delim = ",", col_names = TRUE)

purrr::map(fang_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Fang_class.csv"), delim = ",", col_names = TRUE)

# Gurp --------------------------------------------------------------------
gurp <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Gurp/GSE150724_series_matrix.txt"), getGPL = F)
gurp <- Biobase::phenoData(gurp)@data

gurp_meta <- gurp |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(donor = donor_ch_1,
         enrichment = facs_enrichment_ch_1,
         disease = diagnosis_ch_1) |>
  dplyr::mutate(enrichment = dplyr::case_when(
    stringr::str_detect(enrichment, "gamma- and epsilon") ~ "GammaEpsilon",
    stringr::str_detect(enrichment, "delta") ~ "Delta"),
         sample = base::paste0("Donor", donor, "_", enrichment),
         donor = base::paste0("Donor", donor),
         name = "Gurp",
         donor_study = base::paste0(name, "_", donor),
         age_years = c(48, 48, 52, 52, 51, 51),
         sex = c("m", "m", "m", "m", "m", "m"),
         disease = dplyr::case_when(disease == "non-diabetic" ~ "nd"),
         islet_center = "Integrated Islet Distribution Program (IIDP)") |>
  dplyr::rename(tissue = source_name_ch_1,
                species = organism_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -taxid_ch_1,
                -data_row_count,
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Gurp.csv"), delim = ",", col_names = TRUE)

purrr::map(gurp_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Gurp_class.csv"), delim = ",", col_names = TRUE)

# HPAP --------------------------------------------------------------------
### WHAT ABOUT THE OTHER SHEETS?
# remove t1d donors
hpap <- readxl::read_excel(here::here("islet_cartography_scrna/data_raw/meta/HPAP/Donor_Summary_176.xlsx")) %>%
  dplyr::filter(
    !clinical_diagnosis %in% c(
      "T1DM",
      "T1D control",
      "Recent T1DM Unsuspected",
      "T1DM or MODY, undetermined",
      "T1DM (recent DKA)"
    )
  )

# get names of available sampls
hpap_4 <- utils::read.delim(here::here("islet_cartography_scrna/scripts/download_alignment/HPAP_patch_23/HPAP_patch_23_filt.wget"),
                     header = FALSE) |>
  dplyr::select(donor = V4, sample = V1) |>
  dplyr::distinct() |>
  dplyr::mutate(name = "HPAP_patch_23")

hpap_3 <- utils::read.delim(here::here("islet_cartography_scrna/scripts/download_alignment/HPAP_10x/HPAP_10x.wget"), header = FALSE) |>
  dplyr::select(donor = V1, sample = V1) |>
  dplyr::mutate(dplyr::across(c("donor", "sample"), ~ base::as.character(.))) |>
  dplyr::mutate(name = "HPAP_10x")

hpap_2 <- utils::read.delim(here::here("islet_cartography_scrna/scripts/download_alignment/HPAP_fluidigm/HPAP_fluidigm.wget"),
                     header = FALSE) |>
  dplyr::select(donor = V4, sample = V1) |>
  dplyr::mutate(sample = base::as.character(sample)) |>
  dplyr::mutate(name = "HPAP_fluidigm")

hpap_1 <- utils::read.delim(here::here("islet_cartography_scrna/scripts/download_alignment/HPAP_patch_22/HPAP_patch_22.wget"),
                     header = FALSE) |>
  dplyr::select(donor = V4, sample = V1) |>
  dplyr::distinct() |>
  dplyr::mutate(name = "HPAP_patch_22")

hpap_comb <- purrr::reduce(base::list(hpap_1, hpap_2, hpap_3, hpap_4), dplyr::full_join)

donors <- hpap_comb$donor |> BiocGenerics::unique()

hpap_meta <- hpap %>%
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(donor = donor_id,
                gender = gender,
                disease = clinical_diagnosis,
                hba1c_percent = hba_1_c,
                ethnicity = race,
                disease_duration = disease_duration) %>%
  dplyr::mutate(
    gender = dplyr::case_when(gender == "Male" ~ "m",
                       gender == "Female" ~ "f"),
    disease = dplyr::case_when(disease == "T2DM" ~ "t2d",
                        disease == "T2D control" ~ "nd",
                        disease == "T2DM Gastric bypass" ~ "t2d_gastric_bypass",
                        disease == "T2DM LADA?" ~ "t2d_lada",
                        disease == "T2DM polycystic ovaries" ~ "t2d_polycystic_ovaries",
                        disease == "T2DM (? prediabetic)" ~ "pre"),
    cause_of_death = snakecase::to_snake_case(cause_of_death),
    ethnicity = snakecase::to_snake_case(ethnicity),
    height_cm = base::as.numeric(height_cm),
    weight_kg = base::as.numeric(weight_kg))|>
  dplyr::mutate(dplyr::across(dplyr::where(base::is.character), ~dplyr::na_if(., "NA"))) |>
  dplyr::full_join(y = hpap_comb, by = "donor") |>
  dplyr::filter(donor %in% donors) |>
  dplyr::filter(!age_years < 20) |>
  dplyr::mutate(sample = base::as.character(sample),
                sample = dplyr::case_when(base::is.na(sample) ~ donor,
                                          .default = base::as.character(sample)),
                donor_study = base::paste0(name, "_", donor)) |>
  dplyr::select(dplyr::where(~ !base::all(base::is.na(.)))) %>%
  dplyr::distinct() %>% 
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/HPAP.csv"), delim = ",", col_names = TRUE)

purrr::map(hpap_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  dplyr::mutate(class = dplyr::case_when(class == "POSIXt" ~ "POSIXct",
                .default = base::as.character(class))) %>%
  dplyr::distinct() %>%
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/HPAP_class.csv"), delim = ",", col_names = TRUE)

# Kang --------------------------------------------------------------------
kang<- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Kang/GSE217837_series_matrix.txt"), getGPL = F)
kang <- Biobase::phenoData(kang)@data |>
  dplyr::filter(!base::grepl("Mouse", title))

kang_meta <- kang |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(sample = title,
                disease = condition_ch_1,
                organism = organism_ch_1,
                tissue_source = source_name_ch_1) |>
  dplyr::mutate(
    sample = stringr::str_remove(sample, "Human Pancreatic Islet "),
    sample = stringr::str_remove(sample, "human Pancreatic Islet "),
    sample = stringr::str_replace(sample, " ", "_"),
    disease = dplyr::case_when(
      disease == "Non-diabetic Normal" ~ "nd"
      ),
    donor = dplyr::case_when(
      sample == "scRNA_1" ~ "HP20240",
      sample == "snRNA_1" ~ "HP20240",
      sample == "scRNA_2" ~ "HP20245",
      sample == "snRNA_2" ~ "HP20245",
      sample == "scRNA_3" ~ "HP21189",
      sample == "snRNA_3" ~ "HP21189",
      TRUE ~ "Unknown"
    ),
    age_years = dplyr::case_when(
      donor == "HP20240" ~ 42,
      donor == "HP20245" ~ 24,
      donor == "HP21189" ~ 26,
      TRUE ~ NA_real_
    ),
    sex = dplyr::case_when(
      donor == "HP20240" ~ "m",
      donor == "HP20245" ~ "m",
      donor == "HP21189" ~ "m",
      TRUE ~ NA_character_
    ),
    bmi = dplyr::case_when(
      donor == "HP20240" ~ 23.5,
      donor == "HP20245" ~ 20.0,
      donor == "HP21189" ~ 26.3,
      TRUE ~ NA_real_
    ),
    hba1c_percent = dplyr::case_when(
      donor == "HP20240" ~ 5.4,
      donor == "HP20245" ~ 5.4,
      donor == "HP21189" ~ 5.4,
      TRUE ~ NA_real_
    ),
    islet_center = dplyr::case_when(
      donor == "HP20240" ~ "Prodo Lab (Prodo Aliso Viejo, CA)",
      donor == "HP20245" ~ "Prodo Lab  (Prodo Aliso Viejo, CA)",
      donor == "HP21189" ~ "Prodo Lab  (Prodo Aliso Viejo, CA)",
      TRUE ~ NA_character_
    ),
    history_of_diabetes = dplyr::case_when(
      donor == "HP20240" ~ "No",
      donor == "HP20245" ~ "No",
      donor == "HP21189" ~ "No",
      TRUE ~ NA_character_
    ),
    cause_of_death = dplyr::case_when(
      donor == "HP20240" ~ "stroke",
      donor == "HP20245" ~ "head_trauma",
      donor == "HP21189" ~ "head_trauma",
      TRUE ~ NA_character_
    ),
    estimated_purity_percent = dplyr::case_when(
      donor == "HP20240" ~ 95,
      donor == "HP20245" ~ 95,
      donor == "HP21189" ~ 90,
      TRUE ~ NA_real_
    ),
    estimated_viability_percent = dplyr::case_when(
      donor == "HP20240" ~ 95,
      donor == "HP20245" ~ 95,
      donor == "HP21189" ~ 95,
      TRUE ~ NA_real_
    ),
    culture_time_hours = dplyr::case_when(
      donor == "HP20240" ~ 120,
      donor == "HP20245" ~ 120,
      donor == "HP21189" ~ 120,
      TRUE ~ NA_real_
    ),
    handpicked_to_purity = dplyr::case_when(
      donor == "HP20240" ~ "Yes",
      donor == "HP20245" ~ "Yes",
      donor == "HP21189" ~ "Yes",
      TRUE ~ NA_character_
    ),
    name = "Kang",
    donor_study = base::paste0(name, "_", donor)) |>
      dplyr::rename(tissue = tissue_ch_1,
                    species = organism,
                    molecule = molecule_ch_1) |>
      dplyr::select(-dplyr::contains("data_processing"),
                    -dplyr::contains("characteristics_"),
                    -dplyr::contains("supplementary"),
                    -taxid_ch_1,
                    -data_row_count,
                    -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Kang.csv"), delim = ",", col_names = TRUE)

purrr::map(kang_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Kang_class.csv"), delim = ",", col_names = TRUE)

# Lawlor ------------------------------------------------------------------
lawlor <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Lawlor/GSE86469_series_matrix.txt"), getGPL = F)
lawlor <- Biobase::phenoData(lawlor)@data
# donor info
l_sub_1 <-  readxl::read_xlsx(here::here("islet_cartography_scrna/data_raw/meta/Lawlor/supplemental_tabel_1.xlsx"), skip = 2) |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(donor = unos_id,
                hba1c_percent = hb_a_1_c,
                disease = status,
                age_years = age,
                ethnicity = race) |>
  dplyr::select(-bulk_islet_rna_seq_sample_id, -type_see_fig_1_a) |>
  tidyr::separate_rows(single_cell_preparation, sep = ", ") |>
  dplyr::mutate(dplyr::across(c("sex", "disease"), ~ base::tolower(.)),
                dplyr::mutate(dplyr::across("hba1c_percent", ~ base::as.numeric(.))),
                ethnicity = dplyr::case_when(ethnicity == "AA" ~ "african american",
                                             ethnicity == "W" ~ "white",
                                             ethnicity == "H" ~ "hispanic")) |> 
  dplyr::distinct()

# info about repeated samples
l_sub_2 <- readxl::read_xlsx(here::here("islet_cartography_scrna/data_raw/meta/Lawlor/supplemental_tabel_2.xlsx"), sheet = "Single_cell_stats_and_IDs") |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x))

#combine sample and donor data <- 
lawlor_sample_donor <- l_sub_2 |> 
  dplyr::rename(sample = sampleid) |> 
  dplyr::mutate(single_cell_preparation = case_when(
    grepl("1st", sample) ~ "C1",
    grepl("2nd", sample) ~ "C2",
    grepl("3rd", sample) ~ "C3",
    grepl("4th", sample) ~ "C4",
    grepl("5th", sample) ~ "C5",
    grepl("6th", sample) ~ "C6",
    grepl("7th", sample) ~ "C7",
    grepl("8th", sample) ~ "C8",
    grepl("9th", sample) ~ "C9",
    grepl("10th", sample) ~ "C10",
    grepl("11th", sample) ~ "C11",
    grepl("12th", sample) ~ "C12",
    grepl("13th", sample) ~ "C13")) |> 
  dplyr::select(sample, single_cell_preparation, repeat_sample) |> 
  dplyr::mutate(repeat_sample = snakecase::to_snake_case(repeat_sample)) |> 
  dplyr::full_join(y = l_sub_1, by = "single_cell_preparation") |> 
  dplyr::distinct()

lawlor_meta <- lawlor |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(sample =title,
         sex = "sex_ch_1",
         age_years = "age_ch_1",
         bmi = "bmi_ch_1",
         ethnicity = "race_ch_1",
         inferred_cell_type = "cell_type_ch_1",
         donor = "islet_unos_id_ch_1",
         disease = "disease_ch_1") |>
  dplyr::mutate(dplyr::across(c("bmi", "age_years"), ~ base::as.numeric(.)),
                ethnicity = base::tolower(ethnicity),
                disease = dplyr::case_when(disease == "Non-Diabetic" ~ "nd",
                                           disease == "Type 2 Diabetic" ~ "t2d"),
                sex = dplyr::case_when(sex == "Female" ~ "f",
                                       sex == "Male" ~ "m")) |> 
  dplyr::full_join(y = lawlor_sample_donor) |>
  dplyr::mutate(ethnicity = snakecase::to_snake_case(ethnicity),
                name = "Lawlor",
                donor_study = paste0(name, "_", donor),
                islet_center = "The Integrated Islet Distribution Program (IIDP)",
                on_diabetes_medication = snakecase::to_snake_case(on_diabetes_medication)) |>
  dplyr::rename(tissue = tissue_ch_1,
                species = organism_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -taxid_ch_1,
                -data_row_count,
                -dplyr::contains("category")) |>
  dplyr::distinct() |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Lawlor.csv"), delim = ",", col_names = TRUE)

purrr::map(lawlor_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Lawlor_class.csv"), delim = ",", col_names = TRUE)

# Mauvais_Jarvis ----------------------------------------------------------
mj <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Mauvais_Jarvis/GSE266291_series_matrix.txt"), getGPL = F)
mj<- Biobase::phenoData(mj)@data

mj_meta <- mj |> dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(donor = title,
             sex = "sex_ch_1",
             disease = "disease_state_ch_1") |>
  dplyr::mutate(donor = stringr::str_remove(donor, "Pancreatic islets, non-diabetic, biological replicate "),
         donor = base::paste0("Islet_r", donor),
         sample = donor,
         sex = dplyr::case_when(sex == "Male" ~ "m",
                                sex == "Female" ~ "f"),
         disease = dplyr::case_when(disease == "non-diabetic" ~ "nd"),
         name = "Mauvais_Jarvis",
         donor_study = paste0(name, "_", donor),
         islet_center = "PRODO 559 Laboratories Inc|the Integrated Islet Distribution Program (IIDP)") |>
  dplyr::rename(tissue = tissue_ch_1,
                species = organism_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -taxid_ch_1,
                -data_row_count,
                -dplyr::contains("category"),
                -source_name_ch_1) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Mauvais_Jarvis.csv"), delim = ",", col_names = TRUE)

purrr::map(mj_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Mauvais_Jarvis_class.csv"), delim = ",", col_names = TRUE)

# Motakis -----------------------------------------------------------------
motakis <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Motakis/GSE221156_series_matrix.txt"), getGPL = F)
motakis <- Biobase::phenoData(motakis)@data

motakis_meta <- motakis |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(sample = title,
         sex = sex_ch_1,
         age_years = age_ch_1,
         bmi = bmi_ch_1,
         "hba1c_percent" = hba_1_c_ch_1,
         ethnicity = race_ch_1,
         cause_of_death = cause_of_death_ch_1,
         disease = `disease_state_ch_1`) |>
  dplyr::mutate(name = "Motakis",
                sample = stringr::str_replace(sample, "/", "_"),
                donor = sample,
                sex = dplyr::case_when(sex == "Female" ~ "f",
                                       sex == "Female/Female" ~ "f|f",
                                       sex == "Male" ~ "m",
                                       sex == "Male/Male" ~ "m|m",
                                       sex == "Female/Male" ~ "f|m",
                                       sex == "Male/Female" ~ "m|f"),
         donor_study = paste0(name, "_", donor),
         dplyr::across(c("age_years", "bmi", "hba1c_percent"), ~ base::as.numeric(.)),
         islet_center = "Integrated Islet Distribution Program (IIDP)",
         disease = dplyr::case_when(
           disease %in% c("Type2Diabetic") ~ "t2d",
           disease %in% c("PreDiabetic") ~ "pre",
           disease %in% c("NonDiabetic") ~ "nd",
           disease == "NonDiabetic/PreDiabetic" ~ "nd|pre",
           disease == "NonDiabetic/Type2Diabetic" ~ "nd|t2d",
           disease == "Type2Diabetic/PreDiabetic" ~ "t2d|pre",
           disease == "PreDiabetic/Type2Diabetic" ~ "pre|t2d")) |>
  dplyr::rename(tissue = tissue_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Motakis.csv"), delim = ",", col_names = TRUE)

purrr::map(motakis_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Motakis_class.csv"), delim = ",", col_names = TRUE)

# Muraro ------------------------------------------------------------------
muraro <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Muraro/GSE85241_series_matrix.txt"))
muraro <- Biobase::phenoData(muraro)@data

muraro_meta <- muraro |>
  dplyr::rename_with(~ to_snake_case(.x)) |>
  dplyr::rename(sample = title) |>
  dplyr::mutate(sample = str_remove(sample, ", live sorted cells,"),
         sample = stringr::str_replace_all(sample, " ", "_"),
         donor= stringr::str_extract(sample, "Donor_D[0-9]{2}"),
         name = "Muraro",
         donor_study = base::paste0(name, "_", donor),
         disease = base::rep("nd", length(rownames(muraro))),
         sex = dplyr::case_when(donor == "Donor_D28" ~ "m",
                         donor == "Donor_D29" ~ "m",
                         donor == "Donor_D30" ~ "m",
                         donor == "Donor_D31" ~ "m"),
         bmi = dplyr::case_when(donor == "Donor_D28" ~ 26,
                         donor == "Donor_D29" ~ 22,
                         donor == "Donor_D30" ~ 26,
                         donor == "Donor_D31" ~ 25),
         age_years = dplyr::case_when(donor == "Donor_D28" ~ 54,
                               donor == "Donor_D29" ~ 23,
                               donor == "Donor_D30" ~ 48,
                               donor == "Donor_D31" ~ 59),
         islet_center = "Leiden University Medical Center") |>
  dplyr::rename(species = organism_ch_1,
                tissue = tissue_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Muraro.csv"), delim = ",", col_names = TRUE)

purrr::map(muraro_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Muraro_class.csv"), delim = ",", col_names = TRUE)

# Segerstolpe -------------------------------------------------------------
segerstolpe <- utils::read.delim(here("islet_cartography_scrna/data_raw/meta/Segerstolpe/E-MTAB-5061.sdrf.txt"))

segerstolpe_meta <- segerstolpe  |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(donor = characteristics_individual,
         sample = source_name,
         sex = characteristics_sex,
         age_years = characteristics_age,
         bmi = characteristics_body_mass_index,
         hba1c_percent = characteristics_clinical_information,
         disease = characteristics_disease,
         islet_center = characteristics_biosource_provider,
         tissue = characteristics_cell_type,
         inferred_cell_type = characteristics_inferred_cell_type) |>
  dplyr::mutate(name = "Segerstolpe",
         donor_study = base::paste0(name, "_", donor),
         hba1c_percent = stringr::str_remove(hba1c_percent, "HbA1c"),
         hba1c_percent = stringr::str_remove(hba1c_percent, "%"),
         dplyr::across(c("age_years", "bmi", "hba1c_percent"), ~ base::as.numeric(.)),
         disease = dplyr::case_when(disease == "normal" ~ "nd",
                                    disease == "type II diabetes mellitus" ~ "t2d")) |>
  dplyr::select(-factor_value_single_cell_identifier,
                -factor_value_single_cell_identifier,
                -dplyr::contains("characteristics")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Segerstolpe.csv"), delim = ",", col_names = TRUE)

purrr::map(segerstolpe_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Segerstolpe_class.csv"), delim = ",", col_names = TRUE)

# Shrestha ----------------------------------------------------------------
shrestha <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Shrestha/GSE183568_series_matrix.txt"), getGPL = F)
shrestha <- Biobase::phenoData(shrestha)@data
shrestha_wget_2 <- read.delim(here("islet_cartography_scrna/scripts/download_alignment/Shrestha/Shrestha.wget"), header = F)

shrestha_meta <- shrestha |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(disease = disease_state_ch_1) |>
  dplyr::mutate(age_years = stringr::str_extract(technical_replicates_ch_1, "[0-9]{2}")) %>%
  dplyr::filter(age_years >14) |>
  dplyr::mutate(donor = dplyr::case_when(age_years == 50 ~ "SAMN08768781",
                           age_years == 59 ~ "SAMN08768783",
                           age_years == 66 ~ "R232",
                           age_years == 39 ~ "DON185"),
         treatment = dplyr::case_when(grepl("FACS sorted-Alpha cells", technical_replicates_ch_1) ~ "facs_sorted_alpha_cells",
                               base::grepl("FACS sorted-Beta cells", technical_replicates_ch_1) ~ "facs_sorted_beta_cells"),
         sample = dplyr::case_when(base::grepl("rep1", technical_replicates_ch_1) ~ base::paste0(donor,"_Rep1"),
                                   base::grepl("rep2", technical_replicates_ch_1) ~ base::paste0(donor,"_Rep2"),
                                   base::grepl("rep3", technical_replicates_ch_1) ~ base::paste0(donor,"_Rep3")),
         sample = dplyr::if_else(!base::is.na(treatment), base::paste0(sample, "_facs"), sample),
         sample_id = title,
         preparation = to_snake_case(str_extract(technical_replicates_ch_1, "prep [0-9]")),
         name = "Shrestha",
         islet_center = dplyr::case_when(donor == "SAMN08768781" ~ "Integrated Islet Distribution Program (IIDP)",
                                         donor == "SAMN08768783" ~ "Integrated Islet Distribution Program (IIDP)",
                                         donor == "R232" ~ "Alberta Diabetes Institute (ADI)",
                                         donor == "DON185" ~ "Organ Procurement Organization (OPO)"),
         cause_of_death = dplyr::case_when(donor == "SAMN08768781" ~ "cerebrovascular|stroke",
                                           donor == "SAMN08768783" ~ "cerebrovascular|stroke",
                                           donor == "R232" ~ "cerebrovascular|stroke",
                                           donor == "DON185" ~ "anoxia|drug_intoxication"),
         donor_study = paste0(name, "_", donor),
         cold_ischaemia_time_hours = dplyr::case_when(donor == "SAMN08768781" ~ 9.9,
                                                      donor == "DON185" ~ 8.5),
         estimated_purity_percent = dplyr::case_when(donor == "SAMN08768781" ~ 95,
                                                     donor == "SAMN08768783" ~ 95,
                                                     donor == "R232" ~ 90,
                                                     donor == "DON185" ~ 95),
         culture_time_hours = dplyr::case_when(donor == "SAMN08768781" ~ 19,
                                               donor == "SAMN08768783" ~ 44,
                                               donor == "R232" ~ 144,
                                               donor == "DON185" ~ 36),
         hand_picked_to_purity = "yes",
         other_functional_measurement = "perifusion") |>
  dplyr::mutate(sex = c("m", "m", rep("f", 15)),
         bmi = c(rep(22.4, 2), rep(32.3, 3), rep(18.5, 3), rep(34.76, 9)),
         hba1c_percent = c(rep(NA, 2), rep(NA, 3), rep(6.1, 3), rep(4.7, 9)),
         disease = dplyr::case_when(disease == "Non-diabetic" ~ "nd"),
         dplyr::across(c("age_years", "bmi", "hba1c_percent"), ~ base::as.numeric(.))) |>
  dplyr::filter(sample %in% base::unique(shrestha_wget_2$V1)) |> # only keep samples we sequenced
  dplyr::rename(species = organism_ch_1,
                tissue = tissue_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Shrestha.csv"), delim = ",", col_names = TRUE)


purrr::map(shrestha_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Shrestha_class.csv"), delim = ",", col_names = TRUE)

# Son ---------------------------------------------------------------------
son <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Son/GSE98887_series_matrix.txt"), getGPL = FALSE)
son <- Biobase::phenoData(son)@data

# Create the data frame of donor data provided in supplementary information
son_info <- son_donor_info <- base::data.frame(
  donor = c("ND1", "ND2", "ND3", "ND4", "T2D1", "T2D2", "T2D3", "T2D4", "T2D5", "T2D6"),
  sex = c("m", "f", "f", "f", "m", "m", "m", "f", "m", "m"),
  disease = c("nd", "nd", "nd", "nd", "t2d", "t2d", "t2d", "t2d","t2d","t2d"),
  age_years = c(25, 51, 53, 36, 59, 58, 51, 42, 48, 59),
  bmi = c(26, 20, 34, 24, 32, 39, 25, 28, 44, 33),
  Hba1c_percent = c(NA, NA, 5.3, 5.0, 9.6, 8.9, 6.9, 6.7, 6.6, 6.6),
  treatment = c(NA, NA, NA, NA, "No Tx", "Metformin", "N/A upon admission", "Insulin", "N/A upon admission", "Metformin"),
  disease_duration_years = c(NA, NA, NA, NA, 2, 8, NA, 6, NA, 2)
)
son_meta <- son |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::mutate(sample = stringr::str_split(title, " ", simplify = TRUE)[,1],
                disease = dplyr::case_when(base::grepl("normal", source_name_ch_1) ~ "nd",
                                           base::grepl("type-2 diabetes", source_name_ch_1) ~ "t2d"),
                name = "Son",
                islet_center = "NIH’s Integrated Islet Distribution Program (IIDP)",
                donor = case_when(base::grepl("N01", sample) ~ "ND1",
                                  base::grepl("N02", sample) ~ "ND2",
                                  base::grepl("N03", sample) ~ "ND3",
                                  base::grepl("N04", sample) ~ "ND4",
                                  base::grepl("P06", sample) ~ "T2D1",
                                  base::grepl("P05", sample) ~ "T2D2",
                                  base::grepl("P02", sample) ~ "T2D3",
                                  base::grepl("P03", sample) ~ "T2D4",
                                  base::grepl("P01", sample) ~ "T2D5",
                                  base::grepl("P04", sample) ~ "T2D6",
                                  base::grepl("P07", sample) ~ "excluded"),
                donor_study = base::paste0(name, "_", donor)) |>
  tibble::remove_rownames() |>
  dplyr::left_join(y = son_donor_info, by = "donor") |>
  dplyr::rename(species = organism_ch_1,
                tissue = tissue_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category"),
                -title) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Son.csv"), delim = ",", col_names = TRUE)

purrr::map(son_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Son_class.csv"), delim = ",", col_names = TRUE)

# Tritschler --------------------------------------------------------------
trits <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Tritschler/GSE198623-GPL20301_series_matrix.txt"), getGPL = F)
trits <- Biobase::phenoData(trits)@data

# add data from pdf
trits_meta <- trits |> dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(sample = title) |>
  dplyr::mutate(donor = c("R229", "R239", "R237", "R245", "R266"),
         age_years = c(22, 24, 61, 63, 74),
         sex = c("f", "f", "m", "m", "f"),
         bmi = c(23.0, 22.0, 19.6, 22.3, 29.2),
         hba1c_percent = c(5.3, 5.5, 5.9, 5.6, 6.0),
         name = "Tritschler",
         islet_center = "IsletCore facility (Edmonton, AB, Canada)",
         donor_study = paste0(name, "_", donor)) |>
  dplyr::rename(species = organism_ch_1,
                tissue = source_tissue_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category"),
                -source_name_ch_1) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Tritschler.csv"), delim = ",", col_names = TRUE)

purrr::map(trits_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Tritschler_class.csv"), delim = ",", col_names = TRUE)

# Wang --------------------------------------------------------------------
wang_1 <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Wang/GSE83139-GPL11154_series_matrix.txt"), getGPL = F)
wang_1 <- Biobase::phenoData(wang_1)@data
wang_2 <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Wang/GSE83139-GPL16791_series_matrix.txt"), getGPL = F)
wang_2 <- Biobase::phenoData(wang_2)@data

df_1 <- wang_1 |>
  dplyr::rename(sample = title,
                inferred_cell_type = "curated-cell-type:ch1",
                disease = "health:ch1",
                seq_platform = "instrument_model")

df_2 <- wang_2 |>
  dplyr::rename(sample = title,
                inferred_cell_type = "curated-cell-type:ch1",
                disease = "health:ch1",
                seq_platform = "instrument_model")

wang_meta <- dplyr::bind_rows(df_1, df_2) |>
  dplyr::rename_with(~ snakecase::to_snake_case(.x)) |>
  dplyr::mutate(sample = str_replace(sample, " ", "_"),
                donor = str_extract(sample, "^[^:_]+")) |>
  dplyr::filter(!donor %in% c("ICRH80", "ICRH76", "ACGI428")) |>
  dplyr::mutate(sex = dplyr::case_when(donor == "AAJF122" ~ "m",
                                       donor == "ABAF490" ~ "f",
                                       donor == "ACAP236" ~ "m",
                                       donor == "HP-15041" ~ "m",
                                       donor == "HP-15085-01T2D" ~ "f"),
                cultured_days = dplyr::case_when(donor == "AAJF122" ~ 6,
                                                 donor == "ABAF490" ~ 4,
                                                 donor == "ACAP236" ~ 2,
                                                 donor == "HP-15041" ~ 4,
                                                 base::grepl("HP-15085-01T2D_fresh", sample) ~ 4,
                                                 base::grepl("HP-15085-01T2D_8dcult", sample) ~ 12),
                age_years = dplyr::case_when(donor == "AAJF122" ~ 52,
                                             donor == "ABAF490" ~ 39,
                                             donor == "ACAP236" ~ 21,
                                             donor == "HP-15041" ~ 57,
                                             donor == "HP-15085-01T2D" ~ 37),
                ethnicity = dplyr::case_when(donor == "AAJF122" ~ "asian",
                                             donor == "ABAF490" ~ "white",
                                             donor == "ACAP236" ~ "white",
                                             donor == "HP-15041" ~ "african_american",
                                             donor == "HP-15085-01T2D" ~ "white"),
                bmi = dplyr::case_when(donor == "AAJF122" ~ 29.1,
                                       donor == "ABAF490" ~ 45.2,
                                       donor == "ACAP236" ~ 39,
                                       donor == "HP-15041" ~ 23.98,
                                       donor == "HP-15085-01T2D" ~ 39.3),
                sample = stringr::str_replace(sample, pattern = "::", "_"),
                sample = stringr::str_replace(sample, " ", "_"),
                name = "Wang",
                seq_platform = snakecase::to_snake_case(seq_platform),
                donor_study = base::paste0(name, "_", donor),

                islet_center = "University of Pennsylvania (National Institutes of Health DK-19525)|Integrated Islet Distribution Program (IIDP; http:/iidp.coh.org/)",
                dplyr::across(c("age_years", "bmi", "cultured_days"), ~ base::as.numeric(.))) |>
  dplyr::rename(species = organism_ch_1,
                tissue = tissue_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category"),
                -source_name_ch_1) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Wang.csv"), delim = ",", col_names = TRUE)

purrr::map(wang_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Wang_class.csv"), delim = ",", col_names = TRUE)

# Wang_Sander -------------------------------------------------------------
wang_sander <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Wang_Sander/GSE200044_series_matrix.txt"), getGPL = F)
wang_sander <- Biobase::phenoData(wang_sander)@data
wang_sander_xl <- readxl::read_xlsx(here::here("islet_cartography_scrna/data_raw/meta/Wang_Sander/supplementary_tables.xlsx"),
                                    sheet = "Supplementary Table 1a", skip = 2, col_names = T) |>
  dplyr::rename_with( ~ snakecase::to_snake_case(.x)) |>
  dplyr::select(id, islet_center = human_islet_resource_center, treatment_history, duration_fo_t2d = duration_of_t_2_d)


wang_sander_meta <- wang_sander |> dplyr::filter(molecule_ch1 == "total RNA") |>
  dplyr::rename_with( ~ snakecase::to_snake_case(.x)) |>
  dplyr::rename(
    donor = title,
    id = sample_id_ch_1,
    library_id = library_id_ch_1,
    description_id = description,
    gender = gender_ch_1,
    age_years = age_ch_1,
    bmi = bmi_ch_1,
    hba1c_percent = hba_1_c_ch_1,
    disease = disease_state_ch_1,
    ethnicity = race_ch_1,
    library_id = library_id_ch_1,
    islet_index = islet_index_ch_1,
    purity_percent = purity_ch_1) |>
  dplyr::mutate(
    sample = donor,
    gender = snakecase::to_lower_camel_case(gender),
    disease = case_when(
      disease == "Pre-T2D" ~ "pre",
      disease == "Non-diabetic" ~ "nd",
      disease == "T2D" ~ "t2d"
    ),
    ethnicity = snakecase::to_snake_case(ethnicity),
    description_id = stringr::str_remove(description_id, "snRNA_10x/"),
    name = "Wang_Sander",
    donor_study = base::paste0(name, "_", donor),
    dplyr::across(c("age_years", "bmi", "hba1c_percent"), ~ base::as.numeric(.))
  ) |>
  dplyr::left_join(y = wang_sander_xl, by = "id") |>
  dplyr::rename(species = organism_ch_1,
                tissue = source_name_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Wang_Sander.csv"), delim = ",", col_names = TRUE)

purrr::map(wang_sander_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Wang_Sander_class.csv"), delim = ",", col_names = TRUE)


# Xin ---------------------------------------------------------------------
xin <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Xin/GSE81608_series_matrix.txt"), getGPL = FALSE)
xin <- Biobase::phenoData(xin)@data

xin_meta <- xin |>  dplyr::rename(donor = "donor id:ch1",
               sample = title,
               gender = "gender:ch1",
               age_years = "age:ch1",
               ethnicity = "ethnicity:ch1",
               inferred_cell_type = "cell subtype:ch1",
               disease = "condition:ch1") |>
  dplyr::rename_with( ~ snakecase::to_snake_case(.x)) |>
  dplyr::mutate(donor = stringr::str_replace_all(donor, " ", ""),
         sample = base::paste0(donor, "_", stringr::str_remove(sample, "Pancreatic islet cell sample ")),
         name = "Xin",
         gender = snakecase::to_lower_camel_case(gender),
         ethnicity = case_when(ethnicity == "AA" ~ "african_american",
                               ethnicity == "C" ~ "caucasian",
                               ethnicity == "AI" ~ "asian_indian",
                               ethnicity == "H" ~ "hispanic"),
         cause_of_death = case_when(donor == "NonT2D1" ~ "head_trauma|motor_vehicle_accident",
                                    donor == "NonT2D2" ~ "head_trauma|motor_vehicle_accident",
                                    donor == "NonT2D3" ~ "head_trauma|motor_vehicle_accident",
                                    donor == "NonT2D4" ~ "stroke",
                                    donor == "NonT2D5" ~ "head_trauma|motor_vehicle_accident",
                                    donor == "NonT2D6" ~ "head_trauma|fall",
                                    donor == "NonT2D7" ~ "head_trauma|motor_vehicle_accident",
                                    donor == "NonT2D8" ~ "stroke",
                                    donor == "NonT2D9" ~ "head_trauma|motor_vehicle_accident",
                                    donor == "NonT2D10" ~ "aneurysm",
                                    donor == "NonT2D11" ~ "anorexia",
                                    donor == "NonT2D12" ~ "stroke",
                                    donor == "T2D1" ~ "stroke",
                                    donor == "T2D2" ~ "anorexia",
                                    donor == "T2D3" ~ "anorexia",
                                    donor == "T2D4" ~ "stroke",
                                    donor == "T2D5" ~ "anorexia",
                                    donor == "T2D6" ~ "stroke"),
        hba1c_percent = dplyr::case_when(
           donor == "NonT2D2" ~ 5.1,
           donor == "NonT2D3" ~ 4.9,
           donor == "NonT2D6" ~ 5.3,
           donor == "NonT2D7" ~ 5.1,
           donor == "NonT2D8" ~ 4.9,
           donor == "NonT2D11" ~ 5.4,
           donor == "NonT2D12" ~ 5.2,
           donor == "T2D1" ~ 7.0,
           donor == "T2D2" ~ 7.3,
           donor == "T2D3" ~ 7.4,
           donor == "T2D4" ~ 6.5,
           donor == "T2D5" ~ 6.6,
           donor == "T2D6" ~ 6.9),
         bmi = case_when(donor == "NonT2D1" ~ 21,
                         donor == "NonT2D2" ~ 19,
                         donor == "NonT2D3" ~ 24.5,
                         donor == "NonT2D4" ~ 24.1,
                         donor == "NonT2D5" ~ 31.8,
                         donor == "NonT2D6" ~ 26.7,
                         donor == "NonT2D7" ~ 23.4,
                         donor == "NonT2D8" ~ 27.3,
                         donor == "NonT2D9" ~ 25.4,
                         donor == "NonT2D10" ~ 31.7,
                         donor == "NonT2D11" ~ 28.7,
                         donor == "NonT2D12" ~ 22.8,
                         donor == "T2D1" ~ 24,
                         donor == "T2D2" ~ 39.6,
                         donor == "T2D3" ~ 29.9,
                         donor == "T2D4" ~ 43.1,
                         donor == "T2D5" ~ 43.7,
                         donor == "T2D6" ~ 24.4),
         weight_lbs = case_when(donor == "NonT2D1" ~ 172,
                                donor == "NonT2D2" ~ 118,
                                donor == "NonT2D3" ~ 157,
                                donor == "NonT2D4" ~ 121,
                                donor == "NonT2D5" ~ 215,
                                donor == "NonT2D6" ~ 176,
                                donor == "NonT2D7" ~ 159,
                                donor == "NonT2D8" ~ 170,
                                donor == "NonT2D9" ~ 160,
                                donor == "NonT2D10" ~ 234,
                                donor == "NonT2D11" ~ 142,
                                donor == "NonT2D12" ~ 159,
                                donor == "T2D1" ~ 172,
                                donor == "T2D2" ~ 251,
                                donor == "T2D3" ~ 164,
                                donor == "T2D4" ~ 251,
                                donor == "T2D5" ~ 247,
                                donor == "T2D6" ~ 161),
         disease = dplyr::case_when(disease == "non-diabetic" ~ "nd"),
         islet_center = "Prodo Laboratories, Inc") |>
  dplyr::mutate(donor_study = base::paste0(name, "_", donor),
                dplyr::across(c("age_years", "bmi", "hba1c_percent"), ~ base::as.numeric(.))) |>
  dplyr::rename(species = organism_ch_1,
                tissue = source_name_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category"),
                -tissue_ch_1) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Xin.csv"), delim = ",", col_names = TRUE)

purrr::map(xin_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Xin_class.csv"), delim = ",", col_names = TRUE)

# Xin_Diabetes ------------------------------------------------------------
xin_d <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Xin_Diabetes/GSE114297_series_matrix.txt"), getGPL = FALSE)
xin_d <- Biobase::phenoData(xin_d)@data
xin_d <- xin_d |> dplyr::rename(donor = `individual:ch1`) |>
  dplyr::mutate(donor = str_replace_all(donor, "-", "")) |>
  dplyr::rename_with(~ to_snake_case(.x))

xin_d_sup <- read_xlsx(here("islet_cartography_scrna/data_raw/meta/Xin_Diabetes/supplementary_table_1.xlsx"))
xin_d_sup <- xin_d_sup[-c(13, 14), ]
xin_d_sup_meta <- xin_d_sup |>
  dplyr::rename_with(~ to_snake_case(.x)) |>
  dplyr::rename(donor = donor_id,
                age_years = age) |>
  dplyr::mutate(donor = str_replace_all(donor, "-", ""),
                sample = donor,
         ethnicity = dplyr::case_when(ethnicity == "AA" ~ "african_american",
                               ethnicity == "C" ~ "caucasian",
                               ethnicity == "AI" ~ "asian_indian",
                               ethnicity == "H" ~ "hispanic"),
         cause_of_death = snakecase::to_snake_case(cause_of_death),
         hba1c_percent = stringr::str_extract(hb_a_1_c, "^[^%]+"),
         hba1c_mmol_mol = stringr::str_extract(hb_a_1_c, "(?<=\\()[0-9]+"),
         height_ft = base::as.numeric(stringr::str_extract(height, "^[^']+")),
         height_inches = base::as.numeric(stringr::str_extract(height, "(?<=')[0-9]+")),
         name = "Xin_Diabetes",
         donor_study = base::paste0(name, "_", donor),
         islet_center = "Prodo Laboratories",
         dplyr::across(c("age_years", "bmi", "hba1c_percent", "hba1c_mmol_mol"), ~ base::as.numeric(.))) |>
  dplyr::select(-hb_a_1_c) |>
  dplyr::left_join(y = xin_d, by = "donor") |>
  dplyr::rename(species = organism_ch_1,
                tissue = source_name_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Xin_Diabetes.csv"), delim = ",", col_names = TRUE)

purrr::map(xin_d_sup_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Xin_Diabetes_class.csv"), delim = ",", col_names = TRUE)

# Zhang -------------------------------------------------------------------
zhang <- GEOquery::getGEO(filename = here("islet_cartography_scrna/data_raw/meta/Zhang/GSE195986_series_matrix.txt"),
                getGPL = F)
zhang <- Biobase::phenoData(zhang)@data

zhang_meta <- zhang |> dplyr::rename_with( ~ to_snake_case(.x)) |>
  dplyr::rename(donor = title,
                treatment = characteristics_ch_1) |>
  dplyr::mutate(
    sample = donor,
    hba1c_percent = c(5.8, 5.6, 5.4, 5.8, 5.6, 5.2, 5.4, 6.7, 5.9, 7.6, 6.7),
    bmi = c(31.9, 24.2, 30.7, 35.1, 28.5, 23.5, 23.2, 32.5, 31.7, 24.8, 25.8),
    age_years = c(45, 24, 62, 41, 55, 41, 43, 58, 63, 52, 53),
    name = "Zhang",
    donor_study = base::paste0(name, "_", donor),
    islet_center = "Prodo Laboratories Inc",
    disease = "nd",
    disease = dplyr::case_when(
      base::grepl("T2D", sample) ~ "t2d",
      .default = base::as.character(disease)
    )
  ) |>
  dplyr::filter(!age_years < 20) |>
  dplyr::rename(species = organism_ch_1,
                molecule = molecule_ch_1) |>
  dplyr::select(-dplyr::contains("data_processing"),
                -dplyr::contains("characteristics_"),
                -dplyr::contains("supplementary"),
                -dplyr::contains("category")) |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Zhang.csv"), delim = ",", col_names = TRUE)

purrr::map(zhang_meta, class) |> BiocGenerics::as.data.frame() |> tidyr::pivot_longer(
  cols = dplyr::everything(),
  names_to = "column",
  values_to = "class") |>
  vroom::vroom_write(here::here("islet_cartography_scrna/data/metadata/Zhang_class.csv"), delim = ",", col_names = TRUE)
