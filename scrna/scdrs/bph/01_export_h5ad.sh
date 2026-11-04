cat > convert_to_h5ad.py << 'EOF'
import scanpy as sc
import pandas as pd
import scipy.io as sio
from scipy.sparse import csr_matrix
import os

res_dir = "."

print("=== Info ===")

print("Reading counts_raw.mtx...")
counts = sio.mmread(os.path.join(res_dir, "counts_raw.mtx"))

print("Reading genes and cells Info...")
genes = pd.read_csv(os.path.join(res_dir, "genes.csv"), header=None).iloc[:,0].values
cells = pd.read_csv(os.path.join(res_dir, "cells.csv"), header=None).iloc[:,0].values

print(f"Info: {counts.shape}")
print(f"genes Info: {len(genes)}, cells Info: {len(cells)}")

print("Reading Info data...")
meta = pd.read_csv(os.path.join(res_dir, "metadata.csv"), index_col=0)

counts = csr_matrix(counts)
if counts.shape[0] == len(genes) and counts.shape[1] == len(cells):
    counts = counts.T
    print(f"Info: {counts.shape}")
elif counts.shape[0] == len(cells) and counts.shape[1] == len(genes):
    print(f"Info: {counts.shape}")
else:
    raise ValueError(f"Dimension mismatch: matrix {counts.shape}, genes {len(genes)}, cells {len(cells)}")

adata = sc.AnnData(
    X=counts.astype('float32'),
    obs=meta,
    var=pd.DataFrame(index=genes)
)
adata.obs_names = cells

output_path = os.path.join(res_dir, "scdata.h5ad")
adata.write(output_path)
print(f"\n Info Generating: {output_path}")
print(f"cells Info: {adata.n_obs:,}")
print(f"genes Info: {adata.n_vars:,}")
EOF

python convert_to_h5ad.py
