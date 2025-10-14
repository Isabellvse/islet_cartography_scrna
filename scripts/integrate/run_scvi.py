# run_scvi.py
import os
import sys
import scanpy as sc
import scvi
import pandas as pd

scvi.settings.seed = 0
scvi.settings.dl_num_workers = 90
scvi.settings.dl_persistent_workers = True
print("Last run with scvi-tools version:", scvi.__version__)

# -------------------------
# Arguments
# -------------------------
n = int(sys.argv[1])
adata_file = sys.argv[2]
model_dir  = sys.argv[3]
file_dir   = sys.argv[4]
key_batch0 = sys.argv[5]  # batch_key
key_batch1 = sys.argv[6]  # categorical covariates (or "None")
key_save   = sys.argv[7]
latent_dims = int(sys.argv[8])
embedding_dims = int(sys.argv[9])

# -------------------------
# Load data
# -------------------------
adata = sc.read_h5ad(adata_file)

# -------------------------
# Setup scVI
# -------------------------
categorical_covariates = [key_batch1] if key_batch1 != "None" else None

scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key=key_batch0,
    categorical_covariate_keys=categorical_covariates
)

# -------------------------
# Train model
# -------------------------
model = scvi.model.SCVI(
    adata,
    n_latent=latent_dims,
    batch_representation="embedding",
    batch_embedding_kwargs={'embedding_dim': embedding_dims}
)
model.train()

# -------------------------
# Save model
# -------------------------
model_path = os.path.join(model_dir, f'scvi_{n}_{key_save}')
model.save(model_path, overwrite=True)

# -------------------------
# Save latent embedding
# -------------------------
latent = model.get_latent_representation()
latent_df = pd.DataFrame(
    latent,
    index=adata.obs_names,
    columns=[f"scvi_{i}" for i in range(latent.shape[1])]
)
df_path = os.path.join(file_dir, f"scvi_{n}_{key_save}.csv")
latent_df.to_csv(df_path)

# -------------------------
# Clean up
# -------------------------
del adata
import gc
gc.collect()
