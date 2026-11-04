cat > top1000_only.py << 'EOF'
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import rankdata
import os, sys

TRAIT      = os.environ["TRAIT"]
GROUP_VALUE = os.environ.get("GROUP_VALUE", "BPH")
OUT_TAG = os.environ.get("OUT_TAG", "bph")
H5AD = os.environ.get("H5AD", "scdata.h5ad")
N_TOP, CHUNK = 1000, 500
OUT_ROOT = "./scDRS_out"
SCORE_FILE = f"{OUT_ROOT}/{TRAIT}.full_score.gz"

print(f"Reading {H5AD}...")
adata = sc.read_h5ad(H5AD)
adata_d = adata[adata.obs['group'] == GROUP_VALUE].copy()
print(f"{GROUP_VALUE} cells {adata_d.n_obs}, genes {adata_d.n_vars}")

if not os.path.exists(SCORE_FILE):
    print(f"Not found {SCORE_FILE}, Info Running 1 Info")
    sys.exit(1)

df_full = pd.read_csv(SCORE_FILE, sep='\t', index_col=0)
score = df_full['norm_score'].reindex(adata_d.obs.index)
mask = score.notna().values
print(f"Info {mask.sum()}/{adata_d.n_obs}")
if mask.sum() < adata_d.n_obs * 0.9:
    print("Warning: Info, Info Info")
    sys.exit(1)

s = score.values[mask]
sr = (rankdata(s) - rankdata(s).mean()) / rankdata(s).std()

print(f"Info genes group cor Computing Info ({adata_d.n_vars} genes)...")
all_genes = np.array(adata_d.var_names)
cor_all = np.full(adata_d.n_vars, np.nan)
X_mask = adata_d[mask, :]

for start in range(0, adata_d.n_vars, CHUNK):
    idx = np.arange(start, min(start + CHUNK, adata_d.n_vars))
    Xc = X_mask[:, idx].X
    if hasattr(Xc, 'toarray'):
        Xc = Xc.toarray()
    Xc = np.asarray(Xc, dtype=np.float32)
    for j, gi in enumerate(idx):
        x = Xc[:, j]
        if np.std(x) == 0:
            continue
        xr = rankdata(x)
        xr = (xr - xr.mean()) / xr.std()
        cor_all[gi] = float(np.dot(xr, sr) / len(sr))
    if (start // CHUNK) % 10 == 0:
        print(f"   ... {start}/{adata_d.n_vars}")

df_cor = pd.DataFrame({'gene': all_genes, 'cor': cor_all}).dropna()
df_cor = df_cor.sort_values('cor', ascending=False).reset_index(drop=True)
df_cor['rank'] = df_cor.index + 1
os.makedirs(OUT_ROOT, exist_ok=True)
df_cor.to_csv(f"{OUT_ROOT}/{TRAIT}.{OUT_TAG}_genomewide_cor.csv", index=False)

top = df_cor.head(N_TOP)
top['gene'].to_csv(f"{OUT_ROOT}/{TRAIT}.{OUT_TAG}_top{N_TOP}_genes.txt", index=False, header=False)

print(f"\n Info genes group cor Info: {OUT_ROOT}/{TRAIT}.{OUT_TAG}_genomewide_cor.csv")
print(f"top-{N_TOP} genes Info: {OUT_ROOT}/{TRAIT}.{OUT_TAG}_top{N_TOP}_genes.txt")
print(f"Info threshold cor = {top['cor'].min():.4f}")
print(f"Info 10 Info: {top['gene'].head(10).tolist()}")
EOF

export GROUP_VALUE OUT_TAG H5AD
for TRAIT in ${TRAITS:-BPH_PSA BPH_PSA_6319}; do
  echo ""
  echo "========== top1000: ${TRAIT} =========="
  TRAIT="${TRAIT}" python top1000_only.py
done

echo ""
echo "All completed! Results are in./scDRS_out/"
ls -lh ./scDRS_out/
