# run_scpoli.py
import os
import sys
import scvi
import scanpy as sc
import pandas as pd
from scarches.models.scpoli import scPoli

scvi.settings.seed = 0
scvi.settings.dl_num_workers = 90
scvi.settings.dl_persistent_workers = True  # keeps workers alive between epochs
print("Last run with scvi-tools version:", scvi.__version__)

n = int(sys.argv[1])
adata_file = sys.argv[2]
model_dir  = sys.argv[3]
file_dir   = sys.argv[4]
key_batch = sys.argv[5].split(",")
key_save   = sys.argv[6]
latent_dims = int(sys.argv[7])
embedding_dims = int(sys.argv[8])

adata = sc.read_h5ad(adata_file)

# Move counts to X as dense array
adata.X = adata.layers['counts'].A

# Train scPoli
model = scPoli(adata=adata, condition_keys=key_batch, latent_dim=latent_dims, embedding_dims=embedding_dims)
model.train()

# Save model (optional)
model_path = os.path.join(model_dir, f'scpoli_{n}_{key_save}')
model.save(model_path, save_anndata=True, overwrite=True)

# Save latent embedding
latent = model.get_latent(model.adata, mean=True)
latent_df = pd.DataFrame(latent, index=adata.obs_names,
                         columns=[f"scpoli_{i}" for i in range(latent.shape[1])])
df_path = os.path.join(file_dir, f"scpoli_{n}_{key_save}.csv")
latent_df.to_csv(df_path)

del adata
import gc
gc.collect()
