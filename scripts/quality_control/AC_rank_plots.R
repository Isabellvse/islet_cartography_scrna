# Description -------------------------------------------------------------
# Generate quality control metrics for each study, and save as a .csv file

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(c(here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/")))
set.seed(1000)

# Load --------------------------------------------------------------------
# Ranks
paths <- base::list.files(path = here::here("islet_cartography_scrna/data/quality_control/first_pass/barcode_ranks"), 
                          pattern = "*_ranks.csv", 
                          full.names = TRUE, 
                          recursive = TRUE) |> 
  purrr::set_names(~ stringr::str_extract(.x, pattern = "[^/]+_ic_\\d+_\\d+_\\d+"))

ranks <- purrr::map(paths, ~vroom::vroom(., delim = ",", col_names = TRUE))

# Threshold
paths <- base::list.files(path = here::here("islet_cartography_scrna/data/quality_control/first_pass/barcode_ranks"), 
                          pattern = "*_threshold.csv", 
                          full.names = TRUE, 
                          recursive = TRUE) |> 
  purrr::set_names(~ stringr::str_extract(.x, pattern = "[^/]+_ic_\\d+_\\d+_\\d+"))

thresholds <- purrr::map(paths, ~vroom::vroom(., delim = ",", col_names = TRUE))

# Not empty
paths <- base::list.files(path = here::here("islet_cartography_scrna/data/quality_control/first_pass/barcode_ranks"), 
                          pattern = "*_not_empty", 
                          full.names = TRUE, 
                          recursive = TRUE) |> 
  purrr::set_names(~ stringr::str_extract(.x, pattern = "[^/]+_ic_\\d+_\\d+_\\d+"))

not_empty <- purrr::map(paths, ~vroom::vroom(., delim = ",", col_names = TRUE))

# Rank plot ---------------------------------------------------------------

pdf(file = here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/barcode_rank_plots.pdf"), height = 2, width = 2)
purrr::imap(ranks, function(rank, name){plot_rank(rank = rank, threshold = thresholds[[name]], title = name, non_empty = not_empty[[name]])}) 
dev.off()

#Test one plot
# pdf("test.pdf", height = 2, width = 2)
# plot_rank(rank = ranks[["Wang_Sander_ic_22_16_16"]], threshold = thresholds[["Wang_Sander_ic_22_16_16"]], title = "Wang_Sander_ic_22_16_16", non_empty = not_empty[["Wang_Sander_ic_22_16_16"]])
# dev.off()


# Rank plot as lines ------------------------------------------------------
study_names <- names(ranks) |> stringr::str_remove("_ic_\\d+_\\d+_\\d+") |> unique()
names(study_names) <- study_names
ranks_2 <- purrr::map(study_names, ~ {ranks[grep(.x, names(ranks))] |> 
    dplyr::bind_rows(.id = "name")})

# Define samples that completely fails
failed_samples <- c("Zhang_ic_25_11_11", "Zhang_ic_25_7_7")

pdf(file = here::here("islet_cartography_scrna/data/quality_control/first_pass/plots/barcode_rank_plots_line.pdf"), height = 2, width = 2)
purrr::imap(ranks_2, function(rank, name){plot_rank_line(rank = rank, title = name, failed = failed_samples)}) 
dev.off()
