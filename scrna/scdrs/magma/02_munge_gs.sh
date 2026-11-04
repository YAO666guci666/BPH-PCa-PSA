mkdir -p ./gs_files

rm -f ./gs_files/*.gs

echo "Info PCa_PSA(main analysis)..."
scdrs munge-gs \
    --out-file ./gs_files/PCa_PSA.gs \
    --zscore-file PCa_PSA.tsv \
    --weight zscore \
    --n-max 1000

echo "Info BPH_PSA(main analysis)..."
scdrs munge-gs \
    --out-file ./gs_files/BPH_PSA.gs \
    --zscore-file BPH_PSA.tsv \
    --weight zscore \
    --n-max 1000

echo "Info PCA_PSA_6319(sensitivity analysis)..."
scdrs munge-gs \
    --out-file ./gs_files/PCA_PSA_6319.gs \
    --zscore-file PCA_PSA_6319.tsv \
    --weight zscore \
    --n-max 1000

echo "Info BPH_PSA_6319(sensitivity analysis)..."
scdrs munge-gs \
    --out-file ./gs_files/BPH_PSA_6319.gs \
    --zscore-file BPH_PSA_6319.tsv \
    --weight zscore \
    --n-max 1000

echo "===== Final outputs ====="
ls -l ./gs_files/
echo "All conversions completed!"
