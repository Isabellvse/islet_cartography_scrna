import pandas as pd
import scanpy as sc
import os


# Define if data used to load raw counts should come from Gene or GeneFull_Ex50pAS path
def define_file_path_from_meta(
    metadata_df: pd.DataFrame
):
    # Load the metadata
    # Assuming header=0
    cell_nuclei = metadata_df["cell_nuclei"].unique()

    if cell_nuclei == "cell":
        path = "Gene"
        print(f"[DIR] Path is: ", path)
        return path

    elif cell_nuclei == "nuclei":
        path = "GeneFull_Ex50pAS"
        print(f"[DIR] Subfolder is: ", path)
        return path

    else:
        print("[X] No match of cell or nuclei!")
        return None


def create_anndata_per_sample(
    sample_id: str,
    base_data_dir: str,
    metadata_path: str,
    matrix_filename: str = "matrix.mtx",
    features_filename: str = "features.tsv",
    barcodes_filename: str = "barcodes_prefixed.tsv",
    metadata_barcode_col: int = 6,  # Column index in metadata_df that contains the barcode
    verbose: bool = True,
):
    if verbose:
        print(f" ============= Processing sample: {sample_id} =============")

        # Get subfolder path from meta data
        meta = pd.read_csv(f"{metadata_path}", sep=",", index_col=metadata_barcode_col)

        print(
            f"[OK] 'Meta data' loaded successfully with shape {meta.shape} (rows, columns)"
        )
        print(meta.iloc[0:2, 0:3])

        # subset meta data to only contain sample_id
        meta_sub = meta[meta["sample"] == sample_id]

        # Define subfolder from metadata
        subfolder = define_file_path_from_meta(metadata_df=meta_sub)

        if subfolder is None:
            if verbose:
                print(
                    f"[X] No valid subfolder path found for sample {sample_id}. Check your meta data or sample id! ... Skipping."
                )
                return None  # Early exit

        # Build path to raw data
        sample_raw_data_path = os.path.join(
            base_data_dir, sample_id, "Solo.out", subfolder, "raw"
        )

        if verbose:
            print(f"[DIR] Path to raw data: {sample_raw_data_path}")

            # Check if data path exists for the sample
            if not os.path.exists(sample_raw_data_path):

                if verbose:
                    print(
                        f"[X] Data path not found for {sample_id}: {sample_raw_data_path}. Skipping this sample."
                    )
                    return None

    # Construct paths for matrix, features, and barcodes files
    mtx_file = os.path.join(sample_raw_data_path, matrix_filename)
    print(f"[DIR] MTX file path is: ", mtx_file)
    features_file = os.path.join(sample_raw_data_path, features_filename)
    print(f"[DIR] Feature file path is: ", features_file)
    barcodes_file = os.path.join(sample_raw_data_path, barcodes_filename)
    print(f"[DIR] Barcode file path is: ", barcodes_file)

    # Load barcodes
    bar = pd.read_csv(barcodes_file, sep="\t", index_col=0, header=None)

    if verbose:
        print(
            f"[OK] 'barcodes' loaded for {sample_id} with shape {bar.shape}. Preview:\n{bar.head()}"
        )

    # Read expression matrix
    if verbose:
        print("[..] Reading expression matrix from .mtx file...")
    X_tmp = sc.read_mtx(mtx_file)
    if verbose:
        print(f"[OK] Raw matrix loaded with shape (genes, cells): {X_tmp.X.shape}")

    # Transpose matrix
    if verbose:
        print("[TRA] Transposing matrix to match AnnData format (cells × genes)...")
    X_matrix = X_tmp.X.T
    if verbose:
        print(f"[OK] Matrix transposed. New shape (cells, genes): {X_matrix.shape}")

    # Load features (genes)
    var_df = pd.read_csv(features_file, sep="\t", index_col=0, header=None)

    if verbose:
        print(
            f"[OK] 'features' (genes) loaded with shape {var_df.shape}. Preview:\n{var_df.head()}"
        )

    # Create AnnData object
    if verbose:
        print("[..] Creating AnnData object...")

    X = sc.AnnData(X_matrix, obs=bar, var=var_df)

    if verbose:
        print(f"[OK] AnnData created with shape: {X.shape} (cells × genes)")

    # Subset meta data object to only contain barcodes found in the 'barcodes' file..."
    if verbose:
        print("[SUB] Subsetting metadata to barcodes found in the current sample...")

    meta_bar_subset = meta.loc[[x for x in bar.index.tolist() if x in meta.index]]

    if meta_bar_subset.empty:
        if verbose:
            print(
                "[!!] Subsetting resulted in an empty metadata DataFrame for this sample. This might indicate no matching barcodes. Skipping."
            )
        return None

    if verbose:
        print(f"[OK] Metadata subsetted to {meta_bar_subset.shape} cells.")
        print("[~] Preview of subsetted metadata:")
        print(meta_bar_subset.iloc[0:2, 0:3])

    # Subset AnnData object to only contain barcodes found in the subsetted metadata
    if verbose:
        print("[SUB] Subsetting AnnData object to cells with matching metadata...")

        X_sub = X[meta_bar_subset.index].copy()

    if verbose:
        print(f"[OK] AnnData has been subset to shape: {X_sub.shape} (cells × genes)")

    # Align meta_bar_subset to X_sub.obs's index order
    meta_bar_subset = meta_bar_subset.loc[X_sub.obs.index]

    if verbose:
        assert set(meta_bar_subset.index).issubset(
            set(X_sub.obs.index)
        ), "Some barcodes in obs_bar are not in X_sub.obs"
        print("[OK] Barcodes in meta data are found in AnnData.")

    if verbose:
        # Print which barcodes are missing
        missing = set(meta_bar_subset.index) - set(X.obs.index)

        if missing:
            print(f"Missing barcodes in X.obs: {missing}")
            return None
        else:
            print("[OK] No missing barcodes in meta data.")

            # Align obs_bar to X_sub.obs's index order
            meta_bar_subset = meta_bar_subset.loc[X_sub.obs.index]

            # Check if the indices are now aligned and in the same order
            assert meta_bar_subset.index.equals(
                X_sub.obs.index
            ), "Indices are not in the same order."
            print("[OK] Indices in meta data are in the same order as the AnnData.")

            print("[..] Adding metadata columns from meta data to AnnData...")

            # Ensure obs_bar is aligned to X_sub.obs
            meta_bar_subset = meta_bar_subset.loc[X_sub.obs.index]

            for col in meta_bar_subset.columns:
                X_sub.obs[col] = meta_bar_subset[col]

                print(f"[OK] Added column: '{col}'")

            # Use assert to confirm metadata was copied exactly
            assert X_sub.obs[meta_bar_subset.columns].equals(
                meta_bar_subset
            ), "[X] Metadata was not added correctly to AnnData!"

            print("[OK] All metadata columns added successfully and verified!")
    return X_sub

def get_sample_ids(
    metadata_df: pd.DataFrame
):
    # Get unique sample IDs from the column (e.g., 'sample')
    sample_ids = metadata_df['sample'].unique().tolist()
    print(f"Found {len(sample_ids)} sample IDs:")
    print(sample_ids)
    return sample_ids

    