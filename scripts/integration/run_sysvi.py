# run_sysvi.py
import os
import sys
import scanpy as sc
import scvi
from scvi.external import SysVI
import pandas as pd
import pickle

scvi.settings.seed = 0
scvi.settings.dl_num_workers = 62
scvi.settings.dl_persistent_workers = True  # keeps workers alive between epochs
print("Last run with scvi-tools version:", scvi.__version__)

n = int(sys.argv[1])
adata_file = sys.argv[2]
model_dir  = sys.argv[3]
file_dir   = sys.argv[4]
key_batch0 = sys.argv[5]
key_batch1 = sys.argv[6]
key_save   = sys.argv[7]
latent_dims = int(sys.argv[8])
embedding_dims = int(sys.argv[9])

adata = sc.read_h5ad(adata_file)

# Setup SysVI
SysVI.setup_anndata(
    adata,
    batch_key=key_batch0,
    categorical_covariate_keys=[key_batch1]
)

# SysVI use normalized counts
model = SysVI(
    adata,
    n_latent=latent_dims,
    embed_categorical_covariates=True,
    embedding_kwargs={'embedding_dim': embedding_dims}
)
model.train()

# Save model (optional)
model_path = os.path.join(model_dir, f'sysvi_model_{n}_{key_save}')
model.save(model_path, overwrite=True)

# Save latent embedding
latent = model.get_latent_representation()
latent_df = pd.DataFrame(latent, index=adata.obs_names,
                         columns=[f"sysvi_{i}" for i in range(latent.shape[1])])
df_path = os.path.join(file_dir, f"sysvi_latent_embd_{n}_{key_save}.csv")
latent_df.to_csv(df_path)

# Save model entires add to adata object
sysvi_meta = {
"_sysvi_uuid": adata.uns.get("_sysvi_uuid", None),
"_sysvi_manager_uuid": adata.uns.get("_sysvi_manager_uuid", None)}

uuid_path=os.path.join(file_dir, f'sysvi_{n}_{key_save}_model_uuid.pkl')

with open(uuid_path, "wb") as f:
    pickle.dump(sysvi_meta, f)

del adata
import gc
gc.collect()