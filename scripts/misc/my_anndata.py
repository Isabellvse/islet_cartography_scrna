import pandas as pd
import scanpy as sc
import os
from pathlib import Path
from collections import defaultdict

def vprint(msg, verbose = True):
    """
    Helper function for verbose
    """
    if verbose:
        print(msg)
        
def create_directories(dir_path, verbose):
    """
    Create a directory if it does not already excist

    Args: 
        dir_path: str 
        verbose: bool = True

    Returns:
        None
    """

    if verbose:
        if not os.path.isdir(dir_path):
              os.makedirs(dir_path)
              print(f'{dir_path} Directory created successfully!')
        else:
              print(f'{dir_path} Directory already exists!')

# Define if data used to load raw counts should come from Gene or
# GeneFull_Ex50pAS path
def define_file_path_from_meta(metadata_df, verbose = False):
    """
    Define if raw counts should come from the Gene or GeneFull_Ex50pAS path
    
    Args:
        metadata_df: pd.DataFrame
    Returns:
    str or None: Returns
        - 'Gene' if all entries in 'cell_nuclei' are 'cell'
        - 'GeneFull_Ex50pAS' if all entries in 'cell_nuclei' are 'nuclei'
        - None if the 'cell_nuclei' column contains mixed or unexpected values.
    """
    # Load the metadata
    # Assuming header=0
    cell_nuclei = metadata_df["cell_nuclei"].unique()

    if cell_nuclei == "cell":
        path = "Gene"
        vprint(f"[DIR] Path is: {path}", verbose)
        return path

    elif cell_nuclei == "nuclei":
        path = "GeneFull_Ex50pAS"
        vprint(f"[DIR] Subfolder is: {path}", verbose)
        return path

    else:
        print("[X] No match of cell or nuclei!")
        return None

def create_anndata_per_sample(sample_id, base_data_dir, metadata_path,
    matrix_filename = "matrix.mtx", features_filename = "features.tsv", barcodes_filename = "barcodes_prefixed.tsv", 
    metadata_barcode_col = 6, verbose = True):
    """
    Create an AnnData object for a single sample, where barcodes in meta data file match barcodes in the count files.

    Args:
        sample_id: str,
        base_data_dir: str,
        metadata_path: str,
        matrix_filename: str = "matrix.mtx",
        features_filename: str = "features.tsv",
        barcodes_filename: str = "barcodes_prefixed.tsv",
        metadata_barcode_col: int = 6, # Column index in metadata_df that contains the barcode
        verbose: bool = True
        
    Returns:
        sc.AnnData or None:
        An AnnData object from the scanpy package containing expression data and metadata,
        or None if the processing failed or data is missing.  
    """

    vprint(f" ============= Processing sample: {sample_id} =============", verbose)

    # Get subfolder path from meta data
    meta = pd.read_csv(
        f"{metadata_path}",
        sep=",",
        index_col=metadata_barcode_col)

    vprint(
        f"[OK] 'Meta data' loaded successfully with shape {
            meta.shape} (rows, columns)", verbose)
    vprint(meta.iloc[0:2, 0:3], verbose)

    # subset meta data to only contain sample_id
    meta_sub = meta[meta["sample"] == sample_id]

    # Define subfolder from metadata
    subfolder = define_file_path_from_meta(metadata_df=meta_sub)

    if subfolder is None:
        print(
            f"[X] No valid subfolder path found for sample {sample_id}. Check your meta data or sample id! ... Skipping.")
        return None  # Early exit

    # Build path to raw data
    sample_raw_data_path = os.path.join(
        base_data_dir, sample_id, "Solo.out", subfolder, "raw"
    )

    vprint(f"[DIR] Path to raw data: {sample_raw_data_path}", verbose)

    # Check if data path exists for the sample
    if not os.path.exists(sample_raw_data_path):
        print(
            f"[X] Data path not found for {sample_id}: {sample_raw_data_path}. Skipping this sample.")
        return None

    # Construct paths for matrix, features, and barcodes files
    mtx_file = os.path.join(sample_raw_data_path, matrix_filename)
    features_file = os.path.join(sample_raw_data_path, features_filename)
    barcodes_file = os.path.join(sample_raw_data_path, barcodes_filename)

    # Load barcodes
    bar = pd.read_csv(barcodes_file, sep="\t", index_col=0, header=None)

    vprint(f"[OK] 'barcodes' loaded for {sample_id} with shape {
            bar.shape}. Preview:\n{
            bar.head()}", verbose)

    # Read expression matrix
    vprint("[..] Reading expression matrix from .mtx file...", verbose)
    X_tmp = sc.read_mtx(mtx_file)
    vprint(f"[OK] Raw matrix loaded with shape (genes, cells): {X_tmp.X.shape}", verbose)

    # Transpose matrix
    vprint("[TRA] Transposing matrix to match AnnData format (cells × genes)...", verbose)
    X_matrix = X_tmp.X.T
    vprint(f"[OK] Matrix transposed. New shape (cells, genes): {X_matrix.shape}", verbose)

    # Load features (genes)
    var_df = pd.read_csv(features_file, sep="\t", index_col=0, header=None)
    var_df.columns = ["gene_name", "feature_type"]

    vprint(f"[OK] 'features' (genes) loaded with shape {var_df.shape}. Preview:\n{var_df.head()}", verbose)

    # Create AnnData object
    vprint("[..] Creating AnnData object...", verbose)

    X = sc.AnnData(X_matrix, obs=bar, var=var_df)

    vprint(f"[OK] AnnData created with shape: {X.shape} (cells × genes)", verbose)

    # Subset meta data object to only contain barcodes found in the 'barcodes'
    # file..."
    vprint("[SUB] Subsetting metadata to barcodes found in the current sample...", verbose)

    meta_bar_subset = meta.loc[[
        x for x in bar.index.tolist() if x in meta.index]]

    if meta_bar_subset.empty:
        print("[!!] Subsetting resulted in an empty metadata DataFrame for this sample. This might indicate no matching barcodes. Skipping.")
        return None

    vprint(f"[OK] Metadata subsetted to {meta_bar_subset.shape} cells.", verbose)

    # Subset AnnData object to only contain barcodes found in the subsetted
    # metadata
    vprint("[SUB] Subsetting AnnData object to cells with matching metadata...", verbose)

    X_sub = X[meta_bar_subset.index].copy()

    vprint(
        f"[OK] AnnData has been subset to shape: {
            X_sub.shape} (cells × genes)", verbose)

    # Align meta_bar_subset to X_sub.obs's index order
    meta_bar_subset = meta_bar_subset.loc[X_sub.obs.index]

    assert set(meta_bar_subset.index).issubset(
        set(X_sub.obs.index)
    ), "Some barcodes in obs_bar are not in X_sub.obs"
    vprint("[OK] Barcodes in meta data are found in AnnData.", verbose)

    # Print which barcodes are missing
    missing = set(meta_bar_subset.index) - set(X.obs.index)

    if missing:
        print(f"Missing barcodes in X.obs: {missing}")
        return None
    else:
        vprint("[OK] No missing barcodes in meta data.", verbose)

        # Align obs_bar to X_sub.obs's index order
        meta_bar_subset = meta_bar_subset.loc[X_sub.obs.index]

        # Check if the indices are now aligned and in the same order
        assert meta_bar_subset.index.equals(
            X_sub.obs.index
        ), "Indices are not in the same order."
        vprint("[OK] Indices in meta data are in the same order as the AnnData.", verbose)

        vprint("[..] Adding metadata columns from meta data to AnnData...", verbose)

        # Ensure obs_bar is aligned to X_sub.obs
        meta_bar_subset = meta_bar_subset.loc[X_sub.obs.index]

        for col in meta_bar_subset.columns:
            X_sub.obs[col] = meta_bar_subset[col]
            vprint(f"[OK] Added column: '{col}'", verbose)

        # Use assert to confirm metadata was copied exactly
        assert X_sub.obs[meta_bar_subset.columns].equals(
            meta_bar_subset
        ), "[X] Metadata was not added correctly to AnnData!"

        vprint("[OK] All metadata columns added successfully and verified!", verbose)

    return X_sub


def get_sample_ids(metadata_df):
    """
    Extracts unique sample IDs from the sample column in meta data

    Args:
        metadata_df (pd.DataFrame): A Pandas dataframe
    """
    #
    sample_ids = metadata_df["sample"].unique().tolist()
    return sample_ids


def build_sample_map(metadata_paths, base_data_dirs):
    """
    Match metadata files to their corresponding data directories
    and return a mapping from sample_id to its base_dir and metadata_path.

    Args:
        metadata_paths (List[str]): paths to metadata CSV files
        base_data_dirs (List[str]): paths to preprocessed data directories

    Returns:
        Dict[str, Dict]: sample_id -> { base_data_dir, metadata_path }
        A dictionary keyed by sample IDs (strings), where each value is a dictionary holding info like base data directory and metadata file path.
    """
    # Extract the second to last element in the base_data_dirs (This should be
    # the name of the study, e.g. Segerstolpe)
    base_name_map = {Path(p).parts[-2]: p for p in base_data_dirs}

    # Map sample ids, data paths and meta data paths

    sample_map = {}

    for metadata_path in metadata_paths:
        # Load meta data
        df = pd.read_csv(metadata_path, sep=",", index_col=6)

        # Get sample ids from meta data
        sample_ids = get_sample_ids(df)
        
        # stem = Get the final path component withot suffic (see help(Path.stem))
        # replace = Remove "_metadata" from the final part of the path components
        # (see help(Path.replace))
        metadata_stem = Path(metadata_path).stem.replace("_metadata", "")  # e.g. 'Wang', 'Kang_nuclei'

        # Directly match study name within the base_name_map, this means that we handle cases like 'Wang', 'Wang_Sander'
        # As they are exact matches between meta data name and base_name_map they
        # are not used below.

        if metadata_stem in base_name_map:
            matched_dir = base_name_map[metadata_stem]

        # Special cases - if there is no direct match (e.g., Kang_nuclei → Kang)
        else:
            # Try prefix matching safely (not substring) - "Does the metadata name start with a known folder name?"
            # Kang_nuclei starts with Kang, so we will get this.
            possible_dirs = [
                path
                for name, path in base_name_map.items()
                if metadata_stem.startswith(name)
            ]
            # If we only find one match we continue
            if len(possible_dirs) == 1:
                matched_dir = possible_dirs[0]
            elif len(possible_dirs) > 1:
                # If multiple matches are found, it’s ambiguous
                print(f"[!] Ambiguous match for {metadata_stem}: {possible_dirs}")
                continue
            else:
                # If no matches are found we skip it
                print(f"[!] No match for {metadata_stem}")
                continue

        # Create a named list with sample ids and our data and meta data patches
        for sid in sample_ids:
            sample_map[sid] = {
                "base_data_dir": matched_dir,
                "metadata_path": metadata_path}
        
    return sample_map

def merge_samples_per_study(sample_map, output_dir, subset='yes', verbose = True):
    """
    sample_map: dict
        sample_id -> { base_data_dir, metadata_path }
    output_dir: str
        Where to save combined h5ad files per study (metadata_stem)
    subset: str = 'yes'
        Whether to subset combined AnnData to only include cells marked as "included" (QC filtering)
    """

    # Group samples by study (using metadata_stem)
    study_samples = defaultdict(list)

    for sample_id, info in sample_map.items():
        metadata_stem = Path(info["metadata_path"]).stem.replace("_metadata", "")
        study_samples[metadata_stem].append(sample_id)

    for study, samples in study_samples.items():
        adata_list = []

        for sample_id in samples:
            info = sample_map[sample_id]

            adata = create_anndata_per_sample(
                sample_id=sample_id,
                base_data_dir=info["base_data_dir"],
                metadata_path=info["metadata_path"],
                verbose=verbose
            )

            if adata is not None:
                adata_list.append(adata)

        if adata_list:
            combined_adata = sc.concat(
                adata_list,
                merge="same",
                index_unique=None
            )

            if subset.lower() == 'yes':
                combined_adata = combined_adata[combined_adata.obs.excluded == "included"]

            save_path = os.path.join(output_dir, f"{study}.h5ad")
            combined_adata.write_h5ad(save_path)

            print(f"[OK] Saved combined AnnData for study '{study}' to {save_path}")