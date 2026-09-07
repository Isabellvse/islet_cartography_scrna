# Description -------------------------------------------------------------
# Plot Ucell scores of canonical marker genes scores

# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
base::source(here::here("islet_cartography_scrna/scripts/misc/marker_genes_functions.R"))
set.seed(1000)

iridescent <- khroma::color("iridescent")

# Load --------------------------------------------------------------------
ucell <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/files/ucell_canonical_genes.csv")) |> 
  dplyr::rename(barcode = "...1") |> 
  dplyr::rename(endothelial_islet = islet_endothelial) |> 
  dplyr::select(-acinar_i, -acinar_s, -acinar_endothelial)

# Permutation test --------------------------------------------------------
# Calculate mean original score
orig_score <- ucell |> 
  dplyr::summarise(dplyr::across(.cols = tidyselect::where(is.numeric), .fns = mean,  .names = "{.col}"), 
                   .by = c("cell_type")) |> 
  tidyr::pivot_longer(tidyselect::where(is.numeric), names_to = "gene_set", values_to = paste0("orig_score"))

# Is this cell type's UCell score for this gene set higher than we'd expect by random chance?
# Permutation test
perm_score <- purrr::map(1:1000, function(i) {
  perm_test <- ucell |>
    dplyr::mutate(cell_type = sample(cell_type)) |>
    dplyr::summarise(
      dplyr::across(tidyselect::where(is.numeric), mean),
      .by = cell_type
    ) |>
    tidyr::pivot_longer(-cell_type, names_to = "gene_set", values_to = "perm_score")
  
  is_higher <- orig_score |>
    dplyr::left_join(perm_test, by = c("cell_type", "gene_set")) |>
    dplyr::transmute(is_higher = as.numeric(orig_score <= perm_score)) # drops all not newly made columns
  
  names(is_higher) <- paste0("is_higher_", i)
  return(is_higher)
}) |>
  dplyr::bind_cols()

# add ID columns 
perm_score <- dplyr::bind_cols(orig_score |> dplyr::select(cell_type, gene_set, orig_score), perm_score)

# The star indicates gene lists, which have an fdr less than 0.05, and a scaled ucell score higher than 0.7, 
# there by disregarding gene set scores that are already low in a specific cell type and therefore might become
# significant during the permutation test due to random change anyway. 

test <- perm_score |>
  dplyr::mutate(
    sum   = rowSums(dplyr::across(tidyselect::starts_with("is_higher"))),
    p_val = (sum + 1) / (1000 + 1),
    fdr = stats::p.adjust(p_val, method = "BH")) |> 
  dplyr::group_by(gene_set) |> 
  dplyr::mutate(scaled = (orig_score - min(orig_score))/ (max(orig_score) - min(orig_score)),
                star = dplyr::if_else(fdr <= 0.05 & scaled >= 0.7, "*", "")) |> 
  dplyr::ungroup() |> 
  dplyr::select(cell_type, gene_set, scaled, sum, p_val, fdr, star) 

pdf(
  file = paste0(here::here("islet_cartography_scrna/data/annotate/plot/scaled_ucell_score_heatmap.pdf")),
  height = 2,
  width = 3)
test |> 
  ggplot2::ggplot(ggplot2::aes(y = cell_type, x = gene_set, fill = scaled)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = star), color = "black", size = 1) +
  ggplot2::scale_fill_gradientn(colors = iridescent(5), limits=c(0,1)) +
  ggplot2::labs(y = "Cell type",
                x = "Gene set",
                fill = "UCell score\n(min-max\ncolumns)") +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

