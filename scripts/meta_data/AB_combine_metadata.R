# Description -------------------------------------------------------------
# Here save all the meta data in a list as a .qs file

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/metadata/"))
set.seed(1000)
# load data ---------------------------------------------------------------
# Get paths
files <- base::list.files(path = here::here("islet_cartography_scrna/data/metadata"), pattern = ".csv$", full.names = T)
files_csv <- files[-base::grep("class", files)]
files_class <- files[base::grep("class", files)]

# add names
base::names(files_csv) <- purrr::map_chr(files_csv, ~ stringr::str_extract(string = .x, pattern = "(?<=/)[^/]+(?=\\.csv$)"))
base::names(files_class) <- purrr::map_chr(files_class, ~ stringr::str_extract(string = .x, pattern = "(?<=/)[^/]+(?=\\.csv$)")) |>
  purrr::map_chr(~ stringr::str_remove(string = .x, pattern = "_class"))

# Match order
files_class <- files_class[match(names(files_class), names(files_csv))]

# check all are equal
base::all.equal(base::names(files_csv), base::names(files_class))

# load all dfs
df_loaded <- purrr::map2(files_csv, files_class, ~load_data_with_classes(data_file = .x, class_file = .y))

# split hpap
split_dfs <- df_loaded[["HPAP"]] |> (\(df) base::split(df, factor(df$name)))()
df_loaded[["HPAP"]] <- NULL
df_loaded <- BiocGenerics::append(df_loaded, split_dfs)

# Save --------------------------------------------------------------------
qs2::qs_save(df_loaded, here::here("islet_cartography_scrna/data/metadata/meta_list_1.qs2"))
