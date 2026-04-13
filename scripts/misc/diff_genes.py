## This document will contain code that is used to find differentially expressed genes between clusters
## Packages
import scanpy as sc

import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from joblib import Parallel, delayed

def sample_cell_pseudobulk(ad, sample_key, cluster_key, n_cells, random_state=42):
    """
    Subsample a fixed number of cells per (sample, cluster) group from an AnnData object.

    Parameters
    ----------
    ad : AnnData
        Input AnnData object.
    sample_key : str
        Column name in ad.obs indicating sample identity.
    cluster_key : str
        Column name in ad.obs indicating cluster identity.
    n_cells : int
        Number of cells to sample per (sample, cluster) group.
    random_state : int or np.random.Generator, optional
        Seed or generator for reproducible random subsampling.

    Returns
    -------
    AnnData
        A new AnnData object containing the subsampled cells.
    """
    selected_indices = []

    # Create random number generator
    rng = np.random.default_rng(random_state)

    for (sample_val, cluster_val), indices in ad.obs.groupby([sample_key, cluster_key]).indices.items():
        if len(indices) < n_cells:
            continue
        elif len(indices) == n_cells:
            selected_indices.extend(indices)
        else:
            selected_indices.extend(
                rng.choice(indices, n_cells, replace=False)
            )

    ad_sub = ad[selected_indices].copy()
    return ad_sub

def get_pseudobulk(ad, sample_key, cluster_key, layer="counts", func="sum"):
    """
    Aggregate an AnnData object into pseudobulk profiles by grouping cells
    according to sample and cluster identifiers, then return the aggregated
    AnnData along with count and metadata tables.

    Parameters
    ----------
    ad : AnnData
        Input single-cell AnnData object.
    sample_key : str
        Column in ad.obs specifying sample identity.
    cluster_key : str
        Column in ad.obs specifying cluster or cell-type identity.
    layer : str, optional
        Name of the layer in the AnnData object to aggregate. Defaults to "counts".
    func : str, optional
        Aggregation function to apply when collapsing cells (e.g., "sum", "mean").
        Defaults to "sum".

    Returns
    -------
    dict
        A dictionary containing:
            "pseudobulk_adata" : AnnData
                Aggregated AnnData object.
            "counts" : pandas.DataFrame
                Pseudobulk count matrix with samples × genes.
            "metadata" : pandas.DataFrame
                Corresponding metadata for each pseudobulk sample.
    """
    
    pseudobulk_adata = sc.get.aggregate(
        adata=ad,
        by=[sample_key, cluster_key],
        layer=layer,
        func=func
    )

    counts_df = pd.DataFrame(
        pseudobulk_adata.layers[func],
        columns=pseudobulk_adata.var_names,
        index=pseudobulk_adata.obs_names
    )

    metadata_df = pseudobulk_adata.obs[[sample_key, cluster_key]].copy()
    metadata_df = metadata_df.set_index(counts_df.index)

    output = {
        "pseudobulk_adata": pseudobulk_adata,
        "counts": counts_df,
        "metadata": metadata_df
    }

    return output

def setup_deseq_object(counts, metadata, design, workers=60):
    """
    Initialize a DeseqDataSet object for differential expression analysis.

    Parameters
    ----------
    counts : pandas.DataFrame
        Count matrix with samples as rows and genes as columns.
    metadata : pandas.DataFrame
        Sample-level metadata aligned to the count matrix index.
    design : str, optional
        DESeq2 design formula specifying the model structure.
        An example: "~ ic_id_platform_adjusted_sample + leiden_1.5".
    workers : int, optional
        Number of CPU cores to use for inference. Defaults to 60.

    Returns
    -------
    DeseqDataSet
        A fully initialized DeseqDataSet object configured for downstream DE analysis.
    """
    inference = DefaultInference(n_cpus=workers)

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        inference=inference
    )

    return dds


def prepare_pseudobulk_deseq_analysis(ad, sample_key, cluster_key, n_cells, design, layer="counts", func="sum", random_state = 42, return_all = False, workers=60, no_subset=False):
    """
    This function aggregates single-cell counts into pseudobulk samples based
    on a given sample and cluster key, and then prepares a DESeq2 dataset object
    for downstream differential expression analysis.

    Parameters
    ----------
    ad : AnnData
        Annotated single-cell data object containing expression and metadata.
    n_cells : int
    Number of cells to sample per (sample, cluster) group.
    sample_key : str
        Column in `ad.obs` indicating sample identities.
    cluster_key : str
        Column in `ad.obs` indicating cluster identities.
    design : str
        Design formula for DESeq2 analysis (e.g., "~ condition").
    layer : str, optional (default="counts")
        Which data layer in `ad` to use for generating pseudobulk counts.
    func : str, optional (default="sum")
        Aggregation function to use when summing single-cell counts (e.g., "sum", "mean").
    return_all : bool (default = False) 
        Wether to return all barcodes used to generate pseudobulks
    workers : int, optional (default=60)
        Number of CPU cores to use for parallel computation.
    no_subset : bool, optional (default=False)
        If True, skip subsetting by sample/cluster and generate pseudobulk on all cells.
    
    Returns
    -------
    DESeqDataSet
        A DESeq2-ready dataset object containing pseudobulk counts and metadata
        ready for differential expression analysis.
    
    Notes
    -----
    - This function relies on `sample_cell_pseudobulk` and `get_pseudobulk` for aggregation.
    - Ensure that `ad` is preprocessed (e.g., filtered and normalized) before generating pseudobulk data.
    - The `design` formula should match the experimental setup for proper DESeq2 modeling.
    """
    # Setup n jobs
    inference = DefaultInference(n_cpus=workers)
    
    # Cells for pseudobulk
        # Decide whether to subset cells
    if no_subset:
        ad_sub = ad  # Use full AnnData object
    else:
        ad_sub = sample_cell_pseudobulk(ad=ad, sample_key=sample_key, cluster_key=cluster_key, n_cells=n_cells, random_state=random_state)


    # Pseudobulk
    pseudo_dict = get_pseudobulk(ad = ad_sub, sample_key = sample_key, cluster_key = cluster_key, layer = layer, func = func)

    # Dds object setup
    dds = setup_deseq_object(counts = pseudo_dict["counts"], metadata = pseudo_dict["metadata"], design = design, workers=workers)

    # Fit dispertion and logfoldchanges
    dds.deseq2()

    if return_all:
        return ad_sub.obs[[sample_key, cluster_key]], dds
    else:    
        return dds

def compute_pct_expressing(ad, cluster_key, layer="counts", workers=1):
    """
    Compute percentage of cells expressing each gene (non-zero counts)
    for every cluster in an AnnData object, in parallel. Returns a single
    DataFrame with genes as rows and clusters as columns.

    Parameters
    ----------
    ad : AnnData
        Input AnnData object.
    cluster_key : str
        Column in ad.obs defining cluster identity.
    layer : str, optional
        Layer used for counts. Defaults to "counts".
    n_jobs : int, optional
        Number of parallel workers. Use -1 for all CPUs.

    Returns
    -------
    pandas.DataFrame
        Gene × cluster matrix of percentage-expressing values.
    """

    clusters = ad.obs[cluster_key].unique()

    def compute_for_cluster(cluster):
        ad_c = ad[ad.obs[cluster_key] == cluster]

        X = ad_c.layers[layer]
        if hasattr(X, "toarray"):
            X = X.toarray()

        pct = (X > 0).sum(axis=0) / X.shape[0] * 100
        return pd.Series(
            np.round(pct, 2),
            index=ad_c.var_names,
            name=f"pct_expr_cluster_{cluster}"
        )

    pct_list = Parallel(n_jobs=workers)(
        delayed(compute_for_cluster)(c) for c in clusters
    )

    return pd.concat(pct_list, axis=1)


def diff_genes_two_clusters(dds_obj, cluster_index, cluster_1, cluster_2, workers):
    """
    This function uses the DESeq2-based `DeseqStats` to compare gene expression
    between two specified clusters in a dataset. It runs a Wald test to determine 
    significantly differentially expressed genes.

    Parameters
    ----------
    dds_obj : object
        A DESeqDataSet-like object containing the expression data and metadata.
    cluster_index : str
        The name of the column in the dataset's metadata that indicates cluster labels.
    cluster_1 : str
        The name or label of the first cluster for comparison.
    cluster_2 : str
        The name or label of the second cluster for comparison.
    workers : int
        Number of CPU cores to use for parallel inference.
    
    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the differential expression results for all genes,
        including statistics such as log2 fold change, p-values, and adjusted p-values.
    
    Notes
    -----
    - This function uses the `DefaultInference` class to perform inference in parallel.
    - Make sure `dds_obj` is preprocessed and normalized appropriately before calling this function.
    """
    inference = DefaultInference(n_cpus=workers)
    
    # Run test
    ds = DeseqStats(
        dds_obj,
        contrast=(cluster_index, cluster_1, cluster_2), 
        inference = inference,
        quiet = True)

    # run wald test
    ds.run_wald_test()
    ds.summary()

    # Get results
    results = ds.results_df.copy()

    return results

def extract_results(comparisons, results_list, pct_genes_df):
    # Merge all results
    all_results = pd.concat(results_list, ignore_index=False)
    
    # Merge percentage of genes expressed
    all_results = all_results.merge(pct_genes_df, left_index=True, right_index=True, how="left")
    
    # Convert comparisons to categorical (optional, for sorting later)
    comparison_order = [f"{c1}_vs_{c2}" for c1, c2 in comparisons]
    all_results['comparison'] = pd.Categorical(all_results['comparison'], categories=comparison_order, ordered=True)
    
    # Extract c1 from the comparison string
    all_results['c1'] = all_results['comparison'].str.split('_vs_').str[0]
    
    # Add helper column for percent expressed in c1
    all_results['pct_expr_c1'] = all_results.apply(
        lambda row: row[f"pct_expr_cluster_{row['c1']}"], axis=1
    )
    
    # Reorder columns
    cols_to_move = ['comparison', 'pct_expr_c1', ]
    other_cols = [c for c in all_results.columns if c not in cols_to_move]
    all_results = all_results[cols_to_move + other_cols]

    return(all_results)
    
