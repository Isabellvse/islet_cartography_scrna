# Description -------------------------------------------------------------
# Generate gene sets to calculate gene relevance scores - using genes selected form the PoPs method
# DOI: 10.1038/s41588-023-01443-6

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

base_path <- here::here("islet_cartography_scrna/data/disease_relevance")
dir.create(base_path, showWarnings = FALSE)
dir.create(paste0(base_path, "/", "files"), showWarnings = FALSE)
dir.create(paste0(base_path, "/", "plots"), showWarnings = FALSE)
dir.create(paste0(base_path, "/", "objects"), showWarnings = FALSE)

# Get this https://pmc.ncbi.nlm.nih.gov/articles/instance/10836580/bin/NIHMS1956161-supplement-Supp_Tables.xlsx
# download.file(url = "https://pmc.ncbi.nlm.nih.gov/articles/instance/10836580/bin/NIHMS1956161-supplement-Supp_Tables.xlsx",
#                destfile = paste0(base_path, "/", "files/pops.xlsx"))


# Load --------------------------------------------------------------------
# High confidence genes
traitxgene <- readxl::read_xlsx(here::here("islet_cartography_scrna/data/scdrs/files/NIHMS1956161-supplement-Supp_Tables.xlsx"), 
                                sheet= "S11.High confidence PoPS genes", skip = 1)

# Full name of traits
full_name <- readxl::read_xlsx(here::here("islet_cartography_scrna/data/scdrs/files/NIHMS1956161-supplement-Supp_Tables.xlsx"), 
                               sheet= "S1.Trait summary", skip = 1) |> 
  dplyr::mutate(Description = snakecase::to_snake_case(Description))


# Preprocess --------------------------------------------------------------
# Start with keep metabolic traits, and height 
traits_to_keep <- full_name |> 
  dplyr::filter(Domain == "Metabolic" | Trait == "Height") |> 
  dplyr::pull(Trait)
  
# Filter traits and make negative scores  = 0
traitxgene_flt <- traitxgene |> 
  dplyr::filter(Trait %in% traits_to_keep) |> 
  dplyr::mutate(`PoPS score` = dplyr::if_else(`PoPS score` < 0, 0, `PoPS score`))

stopifnot(traitxgene_flt$`PoPS score` >= 0)
#traitxgene_flt |> ggplot2::ggplot(ggplot2::aes(y = Trait, x = `PoPS score`)) + ggplot2::geom_point()

# with weights

weight <- traitxgene_flt |> dplyr::select(Trait, Gene, score = "PoPS score") |> 
  dplyr::mutate(Gene = base::paste0(Gene, ":", score)) |> 
  dplyr::group_by(Trait) |> 
  dplyr::summarise(GENESET = stringr::str_c(Gene, collapse = ","), .groups = "drop") |> 
  dplyr::left_join(y = full_name) |> 
  dplyr::select(TRAIT = Description, GENESET)

# no weight
noweight <- traitxgene_flt |> 
  dplyr::mutate(`PoPS score` = 1) |> 
  dplyr::select(Trait, Gene, score = "PoPS score") |> 
  dplyr::mutate(Gene = base::paste0(Gene, ":", score)) |> 
  dplyr::group_by(Trait) |> 
  dplyr::summarise(GENESET = stringr::str_c(Gene, collapse = ","), .groups = "drop") |> 
  dplyr::left_join(y = full_name) |> 
  dplyr::select(TRAIT = Description, GENESET) |> 
  dplyr::mutate(TRAIT = paste0("nw_", TRAIT))

# combine
dplyr::bind_rows(weight, noweight) |> 
  vroom::vroom_write(paste0(base_path, "/files/pops_weight_noweight.gs"))

# dplyr::filter(Domain %in% c("Behavioral", "Neurological") | Description %in% c("age_at_menarche", 
#                                                                                  "age_at_menopause", 
#                                                                                  "balding_type_4",
#                                                                                  "loss_of_y",
#                                                                                  "smoking_cigarettes_per_day",
#                                                                                  "smoking_ever_vs_never"))