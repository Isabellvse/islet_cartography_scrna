# Description -------------------------------------------------------------
# Generate quality control metrics for each study, and save as a .csv file

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(c(here::here("islet_cartography_scrna/data/quality_control/first_pass/"),
                   here::here("islet_cartography_scrna/data/quality_control/first_pass/barcode_ranks/"),
                   here::here("islet_cartography_scrna/data/quality_control/first_pass/quality_metrics/")))
set.seed(1000)

# Load --------------------------------------------------------------------
## star_quality control ----
star_quality <- qs2::qs_read(here::here("islet_cartography_scrna/data/quality_control/star_quality_raw.qs2")) |> 
  # Add whether the method is droplet or plate based, ensure library prep is lower case
  purrr::modify_depth(1, ~ dplyr::mutate(., library_prep = base::tolower(library_prep), 
                                         platform = dplyr::case_when(library_prep %in% droplet_based ~ "droplet",
                                                                        library_prep %in% plate_based ~ "plate",
                                                                     library_prep %in% plate_based_bc ~ "plate_barcode")) %>% 
                        dplyr::relocate(platform, .after = "library_prep")) 

# Identify empty droplets -------------------------------------------------
# Donor specific settings
donor_settings_droplet <- list(
  "ic_6_7_7" = list(psi.max = 10), 
  "ic_24_9_9" = list(psi.max = 10),
  "ic_15_26_26" = list(manual_thres = 200),
  "ic_15_27_27" = list(manual_thres = 200),
  "ic_15_28_28" = list(manual_thres = 100),
  "ic_22_13_13" = list(manual_thres = 30),
  "ic_22_16_16" = list(manual_thres = 30),
  "ic_22_9_9" = list(manual_thres = 30),
  "ic_22_7_7" = list(manual_thres = 30),
  "ic_22_8_8" = list(manual_thres = 30),
  "ic_22_5_5" = list(manual_thres = 30)
  )

# out "ic_8_22_22" = list(manual_thres = 1000)
# Get paths for studies which need empty droplet removal
droplet_studies <- star_quality |> 
  purrr::map(function(df){df |> dplyr::filter(platform == "droplet") |> 
      dplyr::pull("name") |> 
      base::unique()}) |>
  purrr::list_c()

# # Get star_quality metrics for to droplet studies
# star_quality_droplet <- star_quality[droplet_studies] 

# Identify empty droplets -------------------------------------------------
# purrr::imap(star_quality_droplet, function(study_metadata, name){
#   base::message("Identifying empty droplets for: ", name)
#   identify_empty_droplets(study_metadata = study_metadata,
#                         study_name = name,
#                         donor_settings = donor_settings_droplet,
#                         save_path = here::here("islet_cartography_scrna/data/quality_control/first_pass/barcode_ranks/"))
#   base::message("Finished ! identifying empty droplets for: ", name)
#   })

# Empty droplets on specific samples --------------------------------------
# Get star_quality metrics for to droplet studies
star_quality_droplet <- star_quality[droplet_studies] |> 
  purrr::modify_depth(1, ~ dplyr::filter(., ic_id %in% "ic_22_16_16"))

i <- star_quality_droplet |> sapply(nrow)>0

star_quality_droplet <- star_quality_droplet[i]

# Identify empty droplets -------------------------------------------------
purrr::imap(star_quality_droplet, function(study_metadata, name){
  base::message("Identifying empty droplets for: ", name)
  identify_empty_droplets(study_metadata = study_metadata,
                          study_name = name,
                          donor_settings = donor_settings_droplet,
                          save_path = here::here("islet_cartography_scrna/data/quality_control/first_pass/barcode_ranks/"))
  base::message("Finished ! identifying empty droplets for: ", name)
})


