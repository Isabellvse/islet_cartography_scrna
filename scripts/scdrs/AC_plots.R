# Description -------------------------------------------------------------
# Plot results from disease relevance

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Load_data ---------------------------------------------------------------
data_disease <- vroom::vroom(here::here("islet_cartography_scrna/data/scdrs/scdrs_metacell_results_merged/group_analysis_merged_cell_type_disease"))
data_disease <- data_disease |> 
  dplyr::mutate(trait2 = factor(trait, levels = c("handedness",
                                                 "hair_colour",
                                                 "waist_cir",
                                                 "bmi",
                                                 "body_fat_percentage",
                                                 "hba1c",
                                                 "t2d",
                                                 "DIAMANTE_T2D",
                                                 "T1D_Forgetta",
                                                 "ins_protein_level",
                                                 "blood_insulin_level",
                                                 "glucagon_measurement",
                                                 "t2d_nephropathy",
                                                 "ischemic_stroke",
                                                 "heart_failure")))

# plot --------------------------------------------------------------------
pdisease <- data_disease |> 
  ggplot2::ggplot(ggplot2::aes(y = trait2, x = group, size = n_fdr_0.1)) +
  ggplot2::geom_point(ggplot2::aes(fill = assoc_mcz), shape = 21) +
  ggplot2::geom_point(data = data_disease[data_disease$assoc_mcp > 0.05,], fill = 'grey', shape = 21) +
  scale_size_continuous(range = c(1, 4), breaks = c(0, 10, 100, 200, 400, 600)) +
  ggplot2::scale_fill_distiller(palette = "RdBu", direction = -1, 
  limits = c(-5, 5),
  oob = scales::squish, 
  name = "MC Z-score\n(grey = p > 0.05)") +
  ggplot2::labs(title = "Disease associations\nat disease and cell type level", x = "Cell type and disease", y = "Trait", 
                size = "Cells with\nFDR > 0.1") +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

pdisease_ht <- data_disease |> 
  ggplot2::ggplot(ggplot2::aes(y = trait2, x = group)) +
  ggplot2::geom_tile(ggplot2::aes(fill = hetero_mcz)) +
  scale_size_continuous(range = c(1, 4), breaks = c(0, 10, 100, 200, 400, 600)) +
  ggplot2::scale_fill_distiller(palette = "RdBu", direction = -1, 
                                limits = c(-5, 5),
                                oob = scales::squish, 
                                name = "MC Z-score\n(grey = p > 0.05)") +
  ggplot2::labs(title = "Disease heterogenity\nat disease and cell type level", x = "Cell type and disease", y = "Trait", 
                size = "Cells with\nFDR > 0.1") +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf(width = 7, height = 3, file = here::here("islet_cartography_scrna/data/scdrs/plots/cell_type_disease_assoc_mcz.pdf"))
pdisease
dev.off()


