traitxgene <- readxl::read_xlsx("C:/Users/isabellvse/Downloads/NIHMS1956161-supplement-Supp_Tables.xlsx", 
                  sheet= "S11.High confidence PoPS genes", skip = 1)

traitxgene |> dplyr::select(TRAIT = Trait, Gene) |> 
  dplyr::group_by(TRAIT) |> 
  dplyr::summarise(GENESET = stringr::str_c(Gene, collapse = ","), .groups = "drop") |> 
  dplyr::mutate(TRAIT = paste0(TRAIT, "_pops")) |> 
  vroom::vroom_write("C:/Users/isabellvse/Downloads/pops_genes.gs")
  
