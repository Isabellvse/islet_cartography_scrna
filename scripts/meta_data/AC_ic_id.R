# Description -------------------------------------------------------------
# Here I generate Islet Cartography (ic) ids
# ic_study_donor_sample
# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/metadata/"))
set.seed(1000)
# load --------------------------------------------------------------------
# study number
study_number <- qs2::qs_read(here::here("islet_cartography_scrna/data/metadata/study_number.qs2"))

# combined meta data frame
meta_list <- qs2::qs_read(here::here("islet_cartography_scrna/data/metadata/meta_list_1.qs2")) |> 
  purrr::modify_depth(1, \(df) dplyr::select(df, name, sample, donor, cell_nuclei, library_prep))

# star quality
# Get paths
files <- base::list.files(path = here::here("islet_cartography_scrna/data/star_quality"), pattern = ".csv$", full.names = T)

# Add names
base::names(files) <- purrr::map_chr(files, ~ stringr::str_extract(string = .x, pattern = "(?<=alignment_qc_)[^/]+(?=\\.csv$)"))

# Load names of samples that have been aligned
star_qc_list <- purrr::imap(files, load_star_quality) |> 
  purrr::modify_depth(1, ~ dplyr::select(.x, name, sample))

# Process -----------------------------------------------------------------
# Get right order of aligned samples
meta_list_order <- meta_list[base::match(base::names(star_qc_list), base::names(meta_list))]
star_qc_list_order <- star_qc_list[base::match(base::names(meta_list), base::names(star_qc_list))]

# Check correct order
base::all.equal(base::names(meta_list_order), base::names(star_qc_list))

# Combine sample donor information with aligned sample information
df_combined_star_list <- purrr::map2(star_qc_list, meta_list_order, ~dplyr::left_join(.x, .y, by = c("sample", "name")))

# check that we don't loose any samples from star 
base::all.equal(base::names(star_qc_list), base::names(df_combined_star_list))
purrr::map2(star_qc_list, df_combined_star_list, ~base::all.equal(base::nrow(.x), base::nrow(.y)))

# is there any na values (any donors and samples that are not in meta)
df_combined_star_list |> 
  purrr::modify_depth(1, ~dplyr::mutate(., dplyr::across(.cols = dplyr::everything(), .fns = is.na)) |> 
                        tidyr::pivot_longer(cols = dplyr::everything(), names_to = "class", values_to = "logical") |> 
                        dplyr::group_by(class, logical) |> 
                        dplyr::tally()) |> 
  dplyr::bind_rows(.id = "study") |> 
  dplyr::filter(logical == "TRUE") 

# 12 samples in Dai does not have any donor information, so they were removed from the meta data, and will also be removed here:
dai_donors_remove <- df_combined_star_list[["Dai"]] |> 
  dplyr::filter(is.na(donor)) |> 
  dplyr::pull("sample")

# remove these cells
df_combined_star_list[["Dai"]] <- df_combined_star_list[["Dai"]] |> 
  dplyr::filter(!sample %in% dai_donors_remove)

# check again if there is there any na values (any donors and samples that are not in meta)
df_combined_star_list |> 
  purrr::modify_depth(1, ~dplyr::mutate(., dplyr::across(.cols = dplyr::everything(), .fns = is.na)) |> 
                        tidyr::pivot_longer(cols = dplyr::everything(), names_to = "class", values_to = "logical") |> 
                        dplyr::group_by(class, logical) |> 
                        dplyr::tally()) |> 
  dplyr::bind_rows(.id = "study") |> 
  dplyr::filter(logical == "TRUE") 

# Combine with overview to get study name and method used

# As the same donors are repeated in the HPAP dataset, we will combine them 
# for generation of ic_ids and then split them again

id_list <- df_combined_star_list |> 
  purrr::modify_depth(1, ~ dplyr::left_join(., y = study_number))

hpap_combined <- dplyr::bind_rows(list(id_list[["HPAP_10x"]], 
                                  id_list[["HPAP_fluidigm"]], 
                                  id_list[["HPAP_patch_22"]],
                                  id_list[["HPAP_patch_23"]])) 
id_list[["HPAP_10x"]] <- NULL
id_list[["HPAP_fluidigm"]] <- NULL
id_list[["HPAP_patch_22"]] <- NULL
id_list[["HPAP_patch_23"]] <- NULL

id_list[["HPAP"]] <- hpap_combined

id_list <- id_list |> purrr::modify_depth(1, ~ dplyr::group_by(., donor) |> 
                                            dplyr::mutate("ic_id_donor" = dplyr::cur_group_id()) %>% 
                                            dplyr::ungroup())

split_dfs <- id_list[["HPAP"]] |> (\(df) base::split(df, factor(df$name)))()
id_list[["HPAP"]] <- NULL
id_list <- BiocGenerics::append(id_list, split_dfs)

id_list <- id_list |> purrr::modify_depth(1, ~ dplyr::group_by(., donor, sample) |> 
                                            dplyr::mutate("ic_id_sample" = dplyr::cur_group_id())  |> 
                                            dplyr::ungroup()  |> 
                                            dplyr::mutate(ic_id = base::paste0("ic", "_", study, "_", ic_id_donor, "_", ic_id_sample)) |> 
                                            dplyr::relocate(ic_id, ic_id_study = study, ic_id_donor, ic_id_sample))
# Test that study ids are correctly assigned

# Extract study ids from study_number
study <- study_number |> dplyr::mutate(pull = base::paste0(name, "_", study)) |> 
                  dplyr::pull("pull")

# Extract study id from id_list
id <- id_list |> purrr::modify_depth(1, ~ dplyr::select(., study = ic_id_study, name)) |> 
  dplyr::bind_rows() |> 
  dplyr::distinct() |> 
  dplyr::mutate(pull = base::paste0(name, "_", study)) |> 
  dplyr::pull("pull")

base::table(study %in% id)

# Save --------------------------------------------------------------------
qs2::qs_save(id_list, here::here("islet_cartography_scrna/data/metadata/id_list.qs2"))
