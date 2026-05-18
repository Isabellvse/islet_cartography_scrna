# Description -------------------------------------------------------------
# Generate gene sets for pops genes

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Get this https://pmc.ncbi.nlm.nih.gov/articles/instance/10836580/bin/NIHMS1956161-supplement-Supp_Tables.xlsx
# download.file(url = "https://pmc.ncbi.nlm.nih.gov/articles/instance/10836580/bin/NIHMS1956161-supplement-Supp_Tables.xlsx",
#               destfile = here::here("islet_cartography_scrna/data/scdrs/files/pops.xlsx"))

traitxgene <- readxl::read_xlsx(here::here("islet_cartography_scrna/data/scdrs/files/NIHMS1956161-supplement-Supp_Tables.xlsx"), 
                                sheet= "S11.High confidence PoPS genes", skip = 1)

full_name <- readxl::read_xlsx(here::here("islet_cartography_scrna/data/scdrs/files/NIHMS1956161-supplement-Supp_Tables.xlsx"), 
                               sheet= "S1.Trait summary", skip = 1) |> 
  dplyr::mutate(Description = snakecase::to_snake_case(Description))

traitxgene |> dplyr::select(Trait, Gene, score = "PoPS score") |> 
  dplyr::mutate(Gene = base::paste0(Gene, ":", score)) |> 
  dplyr::group_by(Trait) |> 
  dplyr::summarise(GENESET = stringr::str_c(Gene, collapse = ","), .groups = "drop") |> 
  dplyr::left_join(y = full_name) |> 
  dplyr::select(TRAIT = Description, GENESET) |> 
  vroom::vroom_write(here::here("islet_cartography_scrna/data/scdrs/gs_files/pops_with_score.gs"))
