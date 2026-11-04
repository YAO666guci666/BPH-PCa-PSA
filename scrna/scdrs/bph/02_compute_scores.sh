TRAITS=("BPH_PSA" "BPH_PSA_6319")
GS_DIR="../magma/gs_files"
OUT_ROOT="./scDRS_out"
mkdir -p "${OUT_ROOT}"

cat > make_unified_bph.py << 'EOF'
import scanpy as sc
import pandas as pd
import os, sys

EXACT_MAP = {
    "b cells": "B Cells", "t cells": "T Cells",
    "basal": "Basal", "club": "Club", "hillock": "Hillock",
    "endothelial": "Endothelial", "malignant": "Malignant",
    "mast cells": "Mast Cells", "myeloid cells": "Myeloid Cells",
    "plasma cells": "Plasma cells",
}

def unify_ct(name):
    key = str(name).strip().lower()
    if 'fibroblast' in key or 'caf' in key:
        return 'Fibroblasts'
    if 'smc' in key or 'pericyte' in key or key == 'smooth muscle':
        return 'Smooth Muscle'
    if 'luminal' in key:
        return 'Luminal'
    return EXACT_MAP.get(key, None)

if os.path.exists("scdata_bph_only.h5ad"):
    try:
        a = sc.read_h5ad("scdata_bph_only.h5ad")
        if "cell_type_unified" in a.obs.columns:
            print("scdata_bph_only.h5ad already exists Info cell_type_unified, skipping Generating")
            print(f"Info cells Info:\n{a.obs['cell_type_unified'].value_counts()}")
            sys.exit(0)
        else:
            print("Warning: Info file Info cell_type_unified columns, rerunning Generating")
    except Exception as e:
        print(f"Warning: Info file Info ({e}), rerunning Generating")

adata = sc.read_h5ad("scdata.h5ad")
for col in ["cell_type_sub", "group"]:
    if col not in adata.obs.columns:
        print(f"obs Info {col} columns, Info columns: {list(adata.obs.columns)}")
        sys.exit(1)

ct = adata.obs["cell_type_sub"].astype(str)
print(f"Info cells Info:\n{ct.value_counts()}\n")

mapped = ct.map(unify_ct)
bad = mapped.isna()
if bad.sum() > 0:
    print("Info cells Info not in Info, Info, Info Running:")
    print(ct[bad].value_counts())
    sys.exit(1)

adata.obs["cell_type_unified"] = pd.Categorical(mapped.values)
print(f"Info cells Info:\n{adata.obs['cell_type_unified'].value_counts()}\n")

adata_bph = adata[adata.obs["group"] == "BPH"].copy()
print(f"BPH group cells Info: {adata_bph.n_obs}")
if adata_bph.n_obs == 0:
    print("BPH group cells Info 0")
    sys.exit(1)

adata_bph.write_h5ad("scdata_bph_only.h5ad")
_ = sc.read_h5ad("scdata_bph_only.h5ad")
print("Info Saved Info validating: scdata_bph_only.h5ad(Info scdata.h5ad Info)")
EOF

echo "Step 1: Generating / validating scdata_bph_only.h5ad..."
python make_unified_bph.py && STEP1_OK=1 || STEP1_OK=0

if [ "${STEP1_OK}" = "1" ]; then
  for TRAIT in "${TRAITS[@]}"; do
    echo ""
    echo "========== Processing trait: ${TRAIT} =========="

    if [ ! -f "${GS_DIR}/${TRAIT}.gs" ]; then
        echo "${GS_DIR}/${TRAIT}.gs does not exist, skipping this trait"
        continue
    fi

    if [ -f "${OUT_ROOT}/${TRAIT}.full_score.gz" ]; then
        echo "${TRAIT}.full_score.gz already exists, skipping scoring"
    else
        echo "Computing ${TRAIT}..."
        scdrs compute-score \
            --h5ad-file ./scdata.h5ad \
            --h5ad-species human \
            --gs-file "${GS_DIR}/${TRAIT}.gs" \
            --gs-species human \
            --out-folder "${OUT_ROOT}/" \
            --n-ctrl 1000 \
            --weight-opt vs

        cd "${OUT_ROOT}"
        [ -f gene_pval.txt ] && mv gene_pval.txt ${TRAIT}_gene_pval.txt
        [ -f gene_score.txt ] && mv gene_score.txt ${TRAIT}_gene_score.txt
        [ -f gene_corr.txt ] && mv gene_corr.txt ${TRAIT}_gene_corr.txt
        if [ -f z_score.full_score.gz ]; then
            mv z_score.full_score.gz ${TRAIT}.full_score.gz
        elif [ -f zscore.full_score.gz ]; then
            mv zscore.full_score.gz ${TRAIT}.full_score.gz
        fi
        if [ -f z_score.score.gz ]; then
            mv z_score.score.gz ${TRAIT}.score.gz
        elif [ -f zscore.score.gz ]; then
            mv zscore.score.gz ${TRAIT}.score.gz
        fi
        cd - > /dev/null
        echo "${TRAIT} Scoring completed"
    fi

    if [ -f "${OUT_ROOT}/${TRAIT}.full_score.gz" ]; then
      echo "${TRAIT} BPH group perform-downstream(unified cell types)..."
      scdrs perform-downstream \
          --h5ad-file ./scdata_bph_only.h5ad \
          --score-file "${OUT_ROOT}/${TRAIT}.full_score.gz" \
          --out-folder "${OUT_ROOT}/" \
          --group-analysis cell_type_unified
      echo "${TRAIT} downstream Completed"
    else
      echo "${TRAIT} Scoring output is missing, skipping downstream"
    fi
  done

  echo ""
  echo "===== All completed, Unified output directory ${OUT_ROOT} ====="
  ls -lh "${OUT_ROOT}"
else
  echo "Step 1 failed, all subsequent steps were skipped. inspect the output above"
fi
