# run_scvi.py
import os
import sys
import scanpy as sc
import pandas as pd
import gc

n = int(sys.argv[1])
adata_file = sys.argv[2]
model_dir  = sys.argv[3]
file_dir   = sys.argv[4]
key_save   = sys.argv[5]
latent_dims = int(sys.argv[6])

adata = sc.read_h5ad(adata_file)

## No integration
sc.tl.pca(adata, n_comps = latent_dims, chunked=True, chunk_size=20000)

# Add X_pca as unintegrated
latent = adata.obsm["X_pca"]
latent_df = pd.DataFrame(latent, index=adata.obs_names,
                         columns=[f"unintegrated_{i}" for i in range(latent.shape[1])])

df_path = os.path.join(file_dir, f"pca_{n}_{key_save}.csv")
latent_df.to_csv(df_path)

del adata
gc.collect()