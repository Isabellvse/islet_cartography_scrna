# Description -------------------------------------------------------------
# The gene length is computed by taking the average length of all transcripts (based on the sum of exon lengths) corresponding to a gene

# Set up ------------------------------------------------------------------
base::source(here::here("islet_cartography_scrna/scripts/misc/set_up.R"))
create_directories(here::here("islet_cartography_scrna/data/quality_control/"))
set.seed(1000)

# Load -------------------------------------------------------------------
gtf_data <- qs2::qs_read(here::here("islet_cartography_scrna/genome_files/gencode.v35.annotation_gtf.qs2"))

# Preprocess -------------------------------------------------------------

# Filter for 'exon' features and calculate the length of each exon
exon_lengths <- gtf_data |>
  dplyr::filter(feature == "exon") |>
  dplyr::mutate(exon_len = end - start + 1)

# Group by transcript and calculate the total length of each transcript
transcript_lengths <- exon_lengths |>
  dplyr::group_by(transcript_id) |>
  dplyr::summarise(transcript_len = base::sum(exon_len))

# Join transcript lengths back to the original data, group by gene_id,
# and calculate the mean transcript length for each gene
gene_lengths <- transcript_lengths |>
  dplyr::left_join(gtf_data |> dplyr::select(transcript_id, gene_id) |> dplyr::distinct(), by = "transcript_id") |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(
    mean_gene_len = base::round(base::mean(transcript_len))
  )

# Save the result to a csv file
vroom::vroom_write(gene_lengths, here::here("islet_cartography_scrna/genome_files/gene_lengths.csv"), delim = ",", col_names = TRUE)
