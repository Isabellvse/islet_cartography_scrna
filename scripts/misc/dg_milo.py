import re
import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData
import pertpy as pt
milo = pt.tl.Milo()
import logging
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

logger = logging.getLogger(__name__)

def milo_da_nhoods(
    mdata: MuData,
    design: str,
    model_contrasts: str | None = None,
    subset_samples: list[str] | None = None,
    add_intercept: bool = True,
    feature_key: str | None = "rna"
) -> AnnData:
    # Extract milo slot
    try:
        sample_adata = mdata["milo"]
    except KeyError:
        logger.error(
            "milo_mdata should be a MuData object with a 'milo' slot - please run milopy.count_nhoods() first"
        )
        raise

    # Extract main feature slot
    adata = mdata[feature_key]

    # Parse covariates from design formula
    covariates = [x.strip() for x in set(re.split(r"\+|\*", design.lstrip("~ ")))]

    # Get sample metadata
    sample_col = sample_adata.uns["sample_col"]
    try:
        sample_obs = adata.obs[covariates + [sample_col]].drop_duplicates()
    except KeyError as e:
        missing_cov = [x for x in covariates if x not in adata.obs.columns]
        logger.warning(f"Covariates {missing_cov} are not columns in adata.obs")
        raise e

    sample_obs.index = sample_obs[sample_col].astype(str)

    # Ensure one-to-one mapping between sample_adata and covariates
    try:
        assert sample_obs.loc[sample_adata.obs_names].shape[0] == len(sample_adata.obs_names)
    except AssertionError:
        logger.warning(
            f"Values in mdata[{feature_key}].obs[{covariates}] cannot be unambiguously assigned to each sample"
        )
        raise

    sample_adata.obs = sample_obs.loc[sample_adata.obs_names]

    # Prepare design dataframe
    design_df = sample_adata.obs[covariates]

    # Get count matrix
    count_mat = sample_adata.X.T.toarray()
    lib_size = count_mat.sum(axis=0)

    # Filter zero-count samples
    keep_smp = lib_size > 0

    # Subset samples if requested
    if subset_samples is not None:
        keep_smp = keep_smp & sample_adata.obs_names.isin(subset_samples)
        design_df = design_df[keep_smp]
        for col in design_df.columns:
            if pd.api.types.is_categorical_dtype(design_df[col]):
                design_df[col] = design_df[col].cat.remove_unused_categories()

    # Filter zero-count neighborhoods
    keep_nhoods = count_mat[:, keep_smp].sum(axis=1) > 0
    counts_filtered = count_mat[np.ix_(keep_nhoods, keep_smp)]
    design_df_filtered = design_df.iloc[keep_smp].copy()

    # Ensure categorical covariates are categories
    for col in design_df_filtered.select_dtypes(exclude=["number"]).columns:
        design_df_filtered[col] = design_df_filtered[col].astype("category")

    design_clean = design if design.startswith("~") else f"~{design}"

    # Build DESeq dataset
    dds = DeseqDataSet(
        counts=pd.DataFrame(counts_filtered.T, index=design_df_filtered.index),
        metadata=design_df_filtered,
        design=design_clean,
        size_factors_fit_type="poscounts",
        refit_cooks=True,
    )
    dds.deseq2()

    # Run contrast or default comparison
    factor_name = design_clean.replace("~", "").split("+")[-1].strip()
    if model_contrasts is not None:
        if "-" not in model_contrasts:
            raise ValueError("Contrast must be in 'GroupA-GroupB' format")
        if "(" in model_contrasts or "+" in model_contrasts.split("-")[1]:
            raise ValueError(
                f"Complex contrasts like '{model_contrasts}' are not supported by pydeseq2"
            )
        group1, group2 = [x.replace(factor_name, "").strip() for x in model_contrasts.split("-")]
        stat_res = DeseqStats(dds, contrast=[factor_name, group1, group2])
    else:
        categories = design_df_filtered[factor_name].cat.categories
        stat_res = DeseqStats(dds, contrast=[factor_name, categories[-1], categories[0]])

    stat_res.summary()
    res = stat_res.results_df.rename(
        columns={"baseMean": "logCPM", "log2FoldChange": "logFC", "pvalue": "PValue", "padj": "FDR"}
    )[['logCPM', 'logFC', 'PValue', 'FDR']]

    # Align results with var
    res.index = sample_adata.var_names[keep_nhoods]
    sample_adata.var = pd.concat([sample_adata.var.drop(res.columns, axis=1, errors='ignore'), res], axis=1)

    # Apply spatial FDR
    milo._graph_spatial_fdr(sample_adata)

    return sample_adata