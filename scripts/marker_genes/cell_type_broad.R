# Description -------------------------------------------------------------
# Identify marker genes from meta analysis of cell types

# Setup -------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
base::source(here::here("islet_cartography_scrna/scripts/misc/marker_genes_functions.R"))
set.seed(1000)
vik <- khroma::color("vik")

library(furrr)
plan(multisession, workers = 64)

base_path <- here::here("islet_cartography_scrna/data/marker_genes")
dir.create(base_path, showWarnings = FALSE)
dir.create(paste0(base_path, "/", "cell_type_broad"), showWarnings = FALSE)
dir.create(paste0(base_path, "/", "cell_type_broad/files"), showWarnings = FALSE)
dir.create(paste0(base_path, "/", "cell_type_broad/plots"), showWarnings = FALSE)


# Load --------------------------------------------------------------------
df_paths <- base::list.files(path = here::here("islet_cartography_scrna/data/annotate/deseq_onevsother_broad_unique"),
                             pattern = ".csv", 
                             full.names = T) |> 
  purrr::set_names(\(vec) (base::basename(vec) |>  
                             stringr::str_remove(".csv")))

df_paths |> 
  purrr::iwalk(\(x, idx) purrr::possibly(process_and_save)(path = x, name = idx, save_path = paste0(base_path, "/cell_type_broad/files")))

# Dot plot ----------------------------------------------------------------
perc <- vroom::vroom(here::here("islet_cartography_scrna/data/annotate/deseq_onevsother_broad_unique/percentage_of_genes_expressed.csv")) |> 
  dplyr::rename_with(~stringr::str_remove(.x, pattern = "pct_expr_cluster_"), 
                     tidyselect::where(is.numeric))

meta_paths <- base::list.files(path = paste0(base_path, "/cell_type_broad/files"),
                               pattern = "*.csv", 
                               full.names = T) |> 
  purrr::set_names(\(vec) (base::basename(vec) |>  
                             stringr::str_remove(".csv")))

# Top 4 markers
markers <- purrr::map(meta_paths, \(df) {
  vroom::vroom(df) |> 
    dplyr::filter(studlab == "metagen (random effects)") |> 
    dplyr::left_join(perc, by = "gene_symbol")
}) |>
  purrr::imap(\(df, name) df |> 
                dplyr::mutate(fdr = stats::p.adjust(pval, method = "BH"))  %>% 
                dplyr::filter(fdr <= 0.01, .data[[name]] >= 50) |> 
                dplyr::select(gene_symbol, TE, seTE, fdr, lower, upper, dplyr::all_of(colnames(perc))) |> 
                dplyr::arrange(desc(TE), desc(.data[[name]])) |>
                head(4)) |> 
  purrr::list_rbind(names_to = "cell_type_broad") |> 
  dplyr::select(gene_symbol, marker = cell_type_broad)

# Percentage long format
perc_long <- perc |> 
  dplyr::filter(gene_symbol %in% markers$gene_symbol) |> 
  tidyr::pivot_longer(-gene_symbol, names_to = "cell_type_broad", values_to = "perc")

logfc <- purrr::map(meta_paths, \(df) {
  vroom::vroom(df) |> 
    dplyr::filter(studlab == "metagen (random effects)") |> 
    dplyr::left_join(perc, by = "gene_symbol")
}) |> 
  purrr::map(\(df) df |> 
               dplyr::filter(studlab == "metagen (random effects)") |> 
               dplyr::filter(gene_symbol %in% markers$gene_symbol) |> 
               dplyr::select(gene_symbol, TE)) |> 
  purrr::list_rbind(names_to = "cell_type_broad") |> 
  dplyr::left_join(markers) |> 
  dplyr::left_join(perc_long)


# dot plot
pdf(
  file = paste0(base_path, "/cell_type_broad/plots/dotplot_top_4_marker_genes.pdf"),
  height = 2,
  width = 4.7)
logfc |>  
  ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder2(gene_symbol, TE, perc), y = cell_type_broad)) +
  ggplot2::geom_point(ggplot2::aes(size = perc, color = TE), shape = 16) +
  ggplot2::scale_size("% Expressed", range = c(0, 2)) +
  ggplot2::scale_color_gradientn(
    colors = vik(256),
    limits = c(-10, 10),
    values = scales::rescale(c(-10, 0, 10)),
    oob = scales::squish, name = "Pooled log2FC"
  ) +
  ggplot2::facet_wrap(~marker, nrow = 1,  scales = "free_x", labeller = labeller(
    marker = function(x)
      ggplot2::label_wrap_gen(width = 10)(
        stringr::str_to_sentence(gsub("_", " ", x))
      ))) +
  ggplot2::scale_y_discrete(
    labels = function(y)(
      stringr::str_to_sentence(gsub("_", " ", y))
    )
  ) +
  my_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                 panel.spacing = unit(0.1, "lines"),
                 axis.title = ggplot2::element_blank())
dev.off()
