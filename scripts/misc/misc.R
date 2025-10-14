
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


# meta data columns -------------------------------------------------------
meta_cols <- c("name",
               "ic_id",
               "ic_id_study",
               "ic_id_donor",
               "donor",
               "sample",
               "age_years",
               "sex",
               "gender",
               "hba1c_percent",
               "Hba1c_percent",
               "bmi",
               "ethnicity",
               "disease",
               "cause_of_death",
               "cell_nuclei",
               "library_prep",
               "islet_center",
               "center_ch_1",
               "organism",
               "cryopreserved",
               "fresh_or_cryo",
               "facs_sorted",
               "facs",
               "sorting_instrument",
               "sequencing_run",
               "cold_ischaemia_time_hours",
               "cold_ischemia_time_h",
               "cold_ischemia_time",
               "cold_ischemia_time_ddhhmmss",
               "culture_time_day",
               "culture_time_hours",
               "cultured_days",
               "geo_accession",
               "cell_type",
               "inferred_cell_type",
               "patched_ch_1", 
               "patched",
               "plate",
               "tissue_source",
               "source_name_ch1",
               "tissue",
               "allocation_via",
               "timefrom_dispersion_days",
               "cultured_days",
               "instrument_model")


# Order of meta data columns ----------------------------------------------
meta_variance_ordered <- c(
  # Study-level meta data and Donor/sample identification
  "ic_id_study", "ic_id_donor", "ic_id_sample", "name", "study", 
  "donor", "sample", "bioproject", "geo_accession", 
  "pmid", "doi", "year_public", "identifier",
  
  # Biometric / clinical data
  "disease", "age_years", "hba_1_c_percent", "sex", "gender", "bmi",
  "ethnicity", "cause_of_death",
  
  # Islet / tissue / culture metadata
  "tissue", "islet_center", 
  "islet_allocation_facility", 
  "islet_culture_medium", "islet_culture_medium_glucose_milimolar", 
  "islet_culture_hours", "cold_ischemia_hours", "islet_fresh_frozen", 
  "islet_isolation_enzyme", "treatment_patch", "treatment_facs",
  
  # Processing / library / sequencing metadata
  "cell_nuclei", "barcode", "study_cell_annotation", "library_prep", 
  "library_layout", "instrument_facs", 
  "instrument_seq", "dissociation_method", "dissociation_tool", 
  "sequencing_run", "strandedness", "type_of_alignment", 
  "star_version",
  
  # Quantification / expression metadata
  "count_quantification", "count_molecule"
)

meta_anndata_ordered <- c(
  # Study-level meta data and Donor/sample identification
  "ic_id_study", "ic_id_donor", "ic_id_sample", "name", "study", 
  "donor", "sample", "bioproject", "geo_accession", 
  "pmid", "doi", "year_public", "identifier",
  
  # Biometric / clinical data
  "disease", "age_years", "hba_1_c_percent", "sex", "gender", "bmi",
  "ethnicity", "ethnicity_broad_harmonized", "ethnicity_sub_harmonized",  
  "cause_of_death", "cause_of_death_broad_harmonized", "cause_of_death_sub_harmonized",
  
  # Islet / tissue / culture metadata
  "tissue", "islet_center", 
  "islet_allocation_facility", 
  "islet_culture_medium", "islet_culture_medium_glucose_milimolar", 
  "islet_culture_hours", "cold_ischemia_hours", "islet_fresh_frozen", 
  "islet_isolation_enzyme", "treatment_patch", "treatment_facs",
  
  # Barcode / cell level meta data 
  "barcode",  "excluded", "n_count", "n_feature",
  "mitochondrial_fraction", "coding_fraction", "contrast_fraction",
  "complexity", "uniquely_mapped_reads_percent",  "unmapped_reads_percent",
  "study_cell_annotation", "study_cell_annotation_harmonized",
  
  # Processing / library / sequencing metadata
  "cell_nuclei", "library_prep", 
  "library_layout", "plate", "instrument_facs", 
  "instrument_seq", "dissociation_method", "dissociation_tool", 
  "sequencing_run", "strandedness", "type_of_alignment", 
  "star_version",
  
  # Quantification / expression metadata
  "rna_count", "count_quantification", "count_molecule"
)