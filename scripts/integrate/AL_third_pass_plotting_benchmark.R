# Description -------------------------------------------------------------
# Here I will plot benchmarking results in a dotplot

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
set.seed(1000)

# Misc --------------------------------------------------------------------
# Define order
metric_type_order = c('Batch correction', 'Bio conservation', 'Aggregate score')
metric_total_order = c('Batch correction', 'Bio conservation', 'Total')

# Define gradient colors
iridescent <- khroma::color("iridescent")

# Load data ---------------------------------------------------------------
per_variable <- vroom::vroom(
  here::here('islet_cartography_scrna/data/integrate/third_pass/files/benchmark_results_per_variable.csv'),
  col_select = c('metric', 'metric_type', 'variable', 'latent', 'score'),
  col_types = vroom::cols(
    metric = vroom::col_character(),
    metric_type = vroom::col_character(),
    variable = vroom::col_character(),
    latent = vroom::col_character(),
    score = vroom::col_double()
  ), 
  .name_repair = snakecase::to_snake_case)

total <- vroom::vroom(
  here::here('islet_cartography_scrna/data/integrate/third_pass/files/benchmark_results_total.csv'),
  col_types = vroom::cols(
    latent = vroom::col_character(),
    metric = vroom::col_character(),
    score = vroom::col_double()
  ),
  .name_repair = snakecase::to_snake_case)

# Define latent order
latent_order <- total |> 
  dplyr::filter(metric == "Total") |> 
  dplyr::arrange(score) |> 
  dplyr::pull(latent)

# Plot --------------------------------------------------------------------
plot_per_variable <- per_variable |> 
  dplyr::mutate(latent = base::factor(latent, levels = latent_order)) |> 
  dplyr::filter(!metric_type %in% c("Batch correction", "Bio conservation")) |> 
  ggplot2::ggplot(ggplot2::aes(x = metric, y = latent)) +
  ggplot2::geom_point(ggplot2::aes(size = score, fill = score), color = "black", shape = 21) +
  ggplot2::scale_fill_gradientn(colors = iridescent(5), limits=c(0,1)) +
  ggplot2::scale_size("Score", range = c(0, 4)) +
  ggplot2::facet_wrap(~variable, nrow = 1, scales = "free_x", labeller = ggplot2:: as_labeller(c("cell_nuclei" = "Cell or Nuclei",
                                                                                                 "ic_id_donor_overall" = "IC ID\ndonor",
                                                                                                 "ic_id_study" = "IC ID\nstudy",
                                                                                                 "library_prep" = "Library\npreparation"))) +
  ggplot2::labs(fill = "Score", size = "Score", title = "Scores per variable") +
  ggplot2::scale_x_discrete(labels = function(x) stringr::word(x, 1)) +
  ggplot2::scale_y_discrete(labels = function(y) stringr::str_replace(y, "X_latent", "IC latent")) +
  my_theme() +
  ggplot2::theme(axis.title = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                 panel.spacing = unit(0.1, "lines"))

plot_total <- total |> 
  dplyr::mutate(latent = base::factor(latent, levels = latent_order),
                metric = base::factor(metric, levels = metric_total_order)) |> 
  ggplot2::ggplot(ggplot2::aes(x = score, y = latent, fill = score)) +
  ggplot2::geom_bar(stat='identity') +
  ggplot2::scale_fill_gradientn(colors = iridescent(5), limits=c(0,1)) +
  ggplot2::facet_wrap(~metric, nrow = 1, labeller = ggplot2:: as_labeller( 
                                                               c("Batch correction" = "Batch",
                                                                 "Bio conservation" = "Bio",
                                                                 "Total" = "Total"))) +
  ggplot2::labs(fill = "Score",
                x = "Score",
                title = "Mean scores\nacross variables") +
  my_theme() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(), 
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                 legend.position = "none",
                 panel.spacing = unit(0.1, "lines")) 

combined_plot <- plot_per_variable + plot_total
combined_plot <- combined_plot +  plot_layout(guides = 'collect', widths = c(0.7, 0.3))

combined_plot

pdf(
  file = here::here("islet_cartography_scrna/data/integrate/third_pass/plot/benchmark_results.pdf"),
  height = 5,
  width = 5
)
combined_plot
dev.off()


