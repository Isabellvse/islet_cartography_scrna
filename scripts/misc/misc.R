# Define droplet and plate based methods ----------------------------------
droplet_based <- c("indrop_v1", "drop_seq", "3p_10x_v3", "3p_10x_v3.1", "3p_10x_v2", "3p_10x_v1", "multiome_10x")
plate_based <- c("smart_seq", "smart_seq2", "smart_seq_ht_ifc", "smarter_seq")
plate_based_bc <- c("sort_seq|cel_seq2") # plate based with barcode

# quality control metrics -------------------------------------------------
qc_metrics <- c("logUMIs", "logFeatures", "mitochondrial_fraction", "ribosomal_fraction", "coding_fraction", "contrast_fraction", "complexity")

