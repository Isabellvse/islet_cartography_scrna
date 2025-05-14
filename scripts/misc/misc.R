# Define droplet and plate based methods ----------------------------------
droplet_based <- c("indrop_v1", "drop_seq", "3p_10x_v3", "3p_10x_v3.1", "3p_10x_v2", "3p_10x_v1", "multiome_10x")
plate_based <- c("smart_seq", "smart_seq2", "smart_seq_ht_ifc", "smarter_seq")
plate_based_bc <- c("sort_seq|cel_seq2") # plate based with barcode

# quality control metrics -------------------------------------------------
qc_metrics_droplet <- c("nUMIs", "nFeatures", "mitochondrial_fraction", "ribosomal_fraction", "coding_fraction", "contrast_fraction", "complexity")
qc_metrics_plate <- c("nCounts", "nFeatures", "mitochondrial_fraction", "ribosomal_fraction", "coding_fraction", "contrast_fraction", "complexity")

qc_met_thres_droplet <- c("nUMIs", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction", "complexity")
qc_met_thres_plate <- c("Uniquely_mapped_reads_%", "Unmapped_reads_%", "nCounts", "nFeatures", "mitochondrial_fraction", "coding_fraction", "contrast_fraction")


qc_star <- c(
  "Number_of_input_reads",
  "Uniquely_mapped_reads_%",
  "Unmapped_reads_%",
  "Average_input_read_length",
  "Average_mapped_length",
  "%_of_chimeric_reads",
  "%_of_reads_mapped_to_multiple_loci",
  "%_of_reads_mapped_to_too_many_loci",
  "%_of_reads_unmapped_other",
  "%_of_reads_unmapped_too_many_mismatches",
  "%_of_reads_unmapped_too_short",
  "Deletion_average_length",
  "Deletion_rate_per_base_%",
  "Insertion_average_length",
  "Insertion_rate_per_base_%",
  "Mismatch_rate_per_base__%",
  "Number_of_chimeric_reads",
  "Number_of_reads_mapped_to_multiple_loci",
  "Number_of_reads_mapped_to_too_many_loci",
  "Number_of_reads_unmapped_other",
  "Number_of_reads_unmapped_too_many_mismatches",
  "Number_of_reads_unmapped_too_short",
  "Number_of_splices_AT/AC",
  "Number_of_splices_Annotated_(sjdb)",
  "Number_of_splices_GC/AG",
  "Number_of_splices_GT/AG",
  "Number_of_splices_Non-canonical",
  "Number_of_splices_Total",
  "Uniquely_mapped_reads_number"
)
