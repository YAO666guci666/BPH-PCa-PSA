cat > calculate_celltype_stats_unified.py << 'EOF'
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
import os, sys

TRAIT        = os.environ["TRAIT"]
GROUP_VALUE  = os.environ.get("GROUP_VALUE", "BPH")
OUT_TAG      = os.environ.get("OUT_TAG", "bph")
CELLTYPE_COL = os.environ.get("CELLTYPE_COL", "cell_type_unified")
H5AD         = os.environ.get("H5AD", "scdata_bph_only.h5ad")
OUT_ROOT     = "./scDRS_out"
SCORE_FILE   = f"{OUT_ROOT}/{TRAIT}.full_score.gz"
GROUP_FILE   = f"{OUT_ROOT}/{TRAIT}.scdrs_group.cell_type_unified"
MC_P         = 0.1

print("Reading data...")
if not os.path.exists(H5AD):
    print(f"Not found {H5AD}, Info Running 1 Info Script")
    sys.exit(1)
adata = sc.read_h5ad(H5AD)

if CELLTYPE_COL not in adata.obs.columns:
    print(f"obs Info {CELLTYPE_COL} columns! Info Running 1 Info Generating unified version h5ad")
    sys.exit(1)

cell_types = adata.obs[CELLTYPE_COL].unique().tolist()
print(f"Info {len(cell_types)} Info cells Info: {cell_types}")
print(f"group columns Info:\n{adata.obs['group'].value_counts()}")

adata_g = adata[adata.obs['group'] == GROUP_VALUE].copy()
print(f"{GROUP_VALUE} group Info cells Info: {adata_g.n_obs}")
if adata_g.n_obs == 0:
    print(f"Info group == '{GROUP_VALUE}' cells, check GROUP_VALUE Info")
    sys.exit(1)

try:
    df_full = pd.read_csv(SCORE_FILE, sep='\t', index_col=0)
    adata_g.obs[f"{TRAIT}_score"] = df_full['norm_score'].reindex(adata_g.obs.index)
    adata_g.obs[f"{TRAIT}_mc_pval"] = df_full['mc_pval'].reindex(adata_g.obs.index)
    n_matched = adata_g.obs[f"{TRAIT}_score"].notna().sum()
    print(f"Info scDRS Info ({n_matched}/{adata_g.n_obs} cells Info)")
    if n_matched < adata_g.n_obs * 0.9:
        print(f"Warning: Info, Info full_score.gz Info cells Info rerunning Running scDRS Generating")
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)

results = []
grouped_scores = []

for cell_type in cell_types:
    mask = adata_g.obs[CELLTYPE_COL] == cell_type
    cells = adata_g[mask]
    if cells.n_obs == 0:
        continue
    scores = cells.obs[f"{TRAIT}_score"].dropna()
    n_cells = len(scores)
    if n_cells == 0:
        continue
    n_sig = (cells.obs[f"{TRAIT}_mc_pval"] < MC_P).sum()
    results.append({
        'cell_type': cell_type,
        'n_cells': n_cells,
        'n_sig': n_sig,
        'prop_sig': n_sig / n_cells,
        'mean_score': scores.mean(),
        'median_score': scores.median(),
        'std_score': scores.std()
    })
    grouped_scores.append(scores.values)

if len(grouped_scores) >= 2:
    kw_stat, kw_pval = stats.kruskal(*grouped_scores)
    print(f"\n Kruskal-Wallis Info: Info cells Info scDRS Info")
    print(f"Info = {kw_stat:.3f}, p = {kw_pval:.4e}")

df_stats = pd.DataFrame(results)

if not os.path.exists(GROUP_FILE):
    print(f"\n Warning: Not found group file: {GROUP_FILE}, skipping group-level FDR")
    df_stats['assoc_mcp'] = np.nan
    df_stats['assoc_fdr'] = np.nan
    df_stats['assoc_mcz'] = np.nan
else:
    gf = pd.read_csv(GROUP_FILE, sep='\t')
    name_col = 'group' if 'group' in gf.columns else gf.columns[0]
    if 'assoc_mcp' not in gf.columns or 'assoc_mcz' not in gf.columns:
        print(f"\n Warning: group file None assoc_mcp or assoc_mcz columns, Info columns: {gf.columns.tolist()}")
        df_stats['assoc_mcp'] = np.nan
        df_stats['assoc_fdr'] = np.nan
        df_stats['assoc_mcz'] = np.nan
    else:
        gf = gf[[name_col, 'assoc_mcp', 'assoc_mcz']].copy()
        gf = gf.sort_values('assoc_mcp').reset_index(drop=True)
        n = len(gf)
        gf['rank_p'] = np.arange(1, n + 1)
        gf['assoc_fdr'] = np.minimum.accumulate((gf['assoc_mcp'] * n / gf['rank_p']).values[::-1])[::-1]
        gf['assoc_fdr'] = np.clip(gf['assoc_fdr'], 0, 1)

        df_stats = df_stats.merge(gf[[name_col, 'assoc_mcp', 'assoc_fdr', 'assoc_mcz']],
                                  left_on='cell_type', right_on=name_col, how='left')
        df_stats = df_stats.drop(columns=[name_col])
        n_na = df_stats['assoc_fdr'].isna().sum()
        print(f"\n group-level FDR and mcz Info Merging ({len(df_stats) - n_na}/{len(df_stats)} cells Info)")
        if n_na > 0:
            print(f"Warning: Info: {df_stats.loc[df_stats['assoc_fdr'].isna(), 'cell_type'].tolist()}")

df_stats = df_stats.sort_values('mean_score', ascending=False)
df_stats['rank'] = range(1, len(df_stats) + 1)

out_cols = ['rank', 'cell_type', 'n_cells', 'n_sig', 'prop_sig',
            'mean_score', 'median_score', 'std_score', 'assoc_mcp', 'assoc_fdr', 'assoc_mcz']
df_stats = df_stats[[c for c in out_cols if c in df_stats.columns]]

os.makedirs(OUT_ROOT, exist_ok=True)
output_file = f"{OUT_ROOT}/{TRAIT}.{OUT_TAG}_cell_type_unified_stats.tsv"
df_stats.to_csv(output_file, sep='\t', index=False)

print(f"\n Saved {GROUP_VALUE} group Info: {output_file}")
print(f"\n {GROUP_VALUE} group Info cells Info scDRS Info (unified cell types):")
display_cols = [c for c in ['rank', 'cell_type', 'n_cells', 'mean_score', 'prop_sig', 'assoc_mcz', 'assoc_fdr'] if c in df_stats.columns]
print(df_stats[display_cols].to_string(index=False))

if 'assoc_fdr' in df_stats.columns and not df_stats['assoc_fdr'].isna().all():
    sig = df_stats[df_stats['assoc_fdr'] < 0.05]
    if len(sig) > 0:
        print(f"\n group-level FDR < 0.05 Info cells Info:")
        for _, r in sig.iterrows():
            print(f"   • {r['cell_type']} (assoc_fdr={r['assoc_fdr']:.4e}, mcz={r['assoc_mcz']:.3f})")
    else:
        print(f"\n Warning: None group-level FDR < 0.05 Info cells Info")

print(f"\n {TRAIT} {GROUP_VALUE} group cells Info Completed")
EOF

export GROUP_VALUE OUT_TAG CELLTYPE_COL H5AD
for TRAIT in ${TRAITS:-BPH_PSA BPH_PSA_6319}; do
  echo ""
  echo "========== Summarizing trait: ${TRAIT} =========="
  TRAIT="${TRAIT}" python calculate_celltype_stats_unified.py
done
