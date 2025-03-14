# Description -------------------------------------------------------------
# Here I generate Islet Cartography (ic) ids
# ic_study_donor_sample
# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/metadata/"))

# load --------------------------------------------------------------------
# data_overview
overview <- readxl::read_xlsx(here::here("islet_cartography_scrna/data/overview_of_data.xlsx"), sheet = "overview", n_max = 25) |> 
  dplyr::select(study, name, library_prep)

# combined meta data frame
meta_list <- qs::qread(here::here("islet_cartography_scrna/data/metadata/meta_list_1.qs")) |> 
  purrr::modify_depth(1, \(df) dplyr::select(df, name, sample, donor))

# star quality
# Get paths
files <- base::list.files(path = here::here("islet_cartography_scrna/data/star_quality"), pattern = ".csv$", full.names = T)

# Add names
base::names(files) <- purrr::map_chr(files, ~ stringr::str_extract(string = .x, pattern = "(?<=alignment_qc_)[^/]+(?=\\.csv$)"))

# Load names of samples that have been aligned
star_qc_list <- purrr::imap(files, load_star_quality) |> 
  purrr::modify_depth(1, ~ dplyr::select(.x, name, sample))

# Process -----------------------------------------------------------------
# Join with overview dataframe, to get method and study number - we will use method later
# df_combined <- df_combined %>%
#   dplyr::left_join(y = overview %>% dplyr::select(study, name, library_prep), by = "name")  |> 
#   dplyr::relocate(name, study)


# Get right order of aligned samples
star_qc_list_order <- star_qc_list[base::match(base::names(meta_list), base::names(star_qc_list))]

# Check correct order
base::all.equal(base::names(star_qc_list_order), base::names(meta_list))

# Combine sample donor information with aligned sample information
df_combined_star_list <- purrr::map2(star_qc_list_order, meta_list, ~dplyr::left_join(.x, .y, by = c("sample", "name")))

# check that we don't loose any samples from star 
base::all.equal(base::names(star_qc_list_order), base::names(df_combined_star_list))
purrr::map2(star_qc_list_order, df_combined_star_list, ~base::all.equal(base::nrow(.x), base::nrow(.y)))

# is there any na values (any donors and samples that are not in meta)
df_combined_star_list |> 
  purrr::modify_depth(1, \(df) dplyr::mutate(df, dplyr::across(.cols = dplyr::everything(), .fns = is.na)) |> 
                        tidyr::pivot_longer(cols = dplyr::everything(), names_to = "class", values_to = "logical") |> 
                        dplyr::group_by(class, logical) |> 
                        dplyr::tally()) |> 
  dplyr::bind_rows(.id = "study") |> 
  dplyr::filter(logical == "TRUE") # No na values, so we should be good to go

# Combine with overview to get study name and method used

# Add IC id
id_list <- df_combined_star_list |> 
  purrr::modify_depth(1, ~ dplyr::left_join(., y = overview) |> 
                        dplyr::group_by(donor) |> 
                        dplyr::mutate("ic_id_donor" = dplyr::cur_group_id()) |> 
                        dplyr::ungroup()  |> 
                        dplyr::group_by(donor, sample) |> 
                        dplyr::mutate("ic_id_sample" = dplyr::cur_group_id())  |> 
                        dplyr::ungroup()  |> 
                        dplyr::mutate(ic_id = base::paste0("ic", "_", study, "_", ic_id_donor, "_", ic_id_sample)) |> 
                        dplyr::relocate(ic_id, ic_id_study = study, ic_id_donor, ic_id_sample))

# Save --------------------------------------------------------------------
qs::qsave(id_list, here::here("islet_cartography_scrna/data/metadata/id_list.qs"))

