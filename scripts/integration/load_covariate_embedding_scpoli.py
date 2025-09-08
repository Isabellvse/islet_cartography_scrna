import os
import sys
import scanpy as sc
import pandas as pd
from pathlib import Path
import subprocess
from scarches.models.scpoli import scPoli

# --- Script part ---
model_path = sys.argv[1]
save_path = sys.argv[2]
embed_name = sys.argv[3]
save_name = sys.argv[4]

model = scPoli.load(model_path)

emb = model.get_conditional_embeddings()[embed_name]

df = pd.DataFrame(
    emb.X.toarray() if hasattr(emb.X, "toarray") else emb.X,
    index=emb.obs_names,
    columns=emb.var_names
)

df_path = os.path.join(save_path, f"{save_name}.csv")
df.to_csv(df_path)