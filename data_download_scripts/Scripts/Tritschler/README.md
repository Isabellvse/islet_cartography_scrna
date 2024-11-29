## Tritschler

- `GEO`: GSE198623
- `PMID`: 36113773
- `Tissue`: Isolated islets


### Dataset at a glance

- `n samples`: 5 Non diabetic
- `method`: Single-Cell RNA-seq, 10x Genomics Single Cell 3' v2

### Alignment paramters

- `--soloType CB_UMI_Simple:` Sets the type of single-cell RNA-seq data, in this case, a simple cell barcode (CB) and unique molecular identifier (UMI) setup.
- `--soloFeatures Gene GeneFull:` Specifies the features to be counted, including genes and full gene bodies.
- `--soloCellFilter None:` Disables cell filtering, meaning all barcodes will be considered.
- `--soloUMIlen 10:` Sets the length of the UMI to 10 bases.
- `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts:` Allows for one mismatch and multiple N bases in the cell barcode matching against the whitelist.
- `--clipAdapterType CellRanger4:` Uses the adapter clipping method similar to Cell Ranger version 4.
- `--outFilterScoreMin 30:` Sets the minimum score for filtering out low-quality reads.
- `--soloMultiMappers EM:` Uses the Expectation-Maximization algorithm for handling multi-mapping reads.
- `--soloUMIfiltering MultiGeneUMI_CR:` Filters UMIs that map to multiple genes, similar to Cell Ranger.
- `--soloUMIdedup 1MM_CR:` Deduplicates UMIs allowing for one mismatch, similar to Cell Ranger.
- `--outMultimapperOrder Random:` Outputs multi-mappers in a random order.
- `--outSAMmultNmax 1:` Limits the number of alignments per read in the output SAM file to 1.
- `--soloCBwhitelist 737K-august-2016.txt:` Specifies the file containing the whitelist of valid cell barcodes.
- `--soloBarcodeReadLength 0:` Uses the full length of the barcode read.