# Code Verification and Revision Record

## Actual workflow, naming, and deduplication

- The main analysis uses three binary GWAS datasets, and continuous PSA is a separate sensitivity analysis. `*_MTAG.txt` was only a legacy filename and has been replaced by the original GWAS aliases.
- All directories and scripts use English names. The main analysis is numbered continuously from `01_preprocess_gwas` to `36_table1_colocalization`, and the continuous-PSA sensitivity analysis from `01_preprocess_gwas` to `16_hdl`.
- Figure 5 depends on LAVA, so the final order is fixed as `31_placo -> 32_prepare_placo_fuma -> 33_lava -> 34_placo_gene_figure -> 35_open4gene -> 36_table1_colocalization`.
- Shared main-analysis and continuous-PSA code is stored only in the main-analysis directories. Short sensitivity entry scripts set `GCST90476319_PSA` and call the corresponding main scripts.
- The three scoring workflows differ only in inputs, thresholds, and output descriptions and share `score_core.R` through three short entry scripts.
- BPH and PCa share the same h5ad conversion, cell-type statistics, and top-1000 code. Step 02 remains separate because the cell sets and downstream settings differ.
- A shared subcluster merge-back script writes `cell_type_sub` by cell barcode, preserves the original UMAP, adds `umap_new`, and outputs `newBPH.rds` or `newPCa.rds`.
- Internal input and output names for scoring, HyPrColoc, Table 1, and logs are in English. Former pinyin names such as `PLACO_JIYIN` have been renamed to `placo_gene_figure`.
- Empty directory placeholders have been removed from the archive, and all documentation and manifest files are in English.

## Corrections implemented

1. Restored the three actual inputs in main-analysis steps 01/02 and removed unrelated disease examples; removed case/control fields for continuous PSA.
2. Removed unrelated polyp inputs and standardized the PCa single-trait filename as `ieu-b-4809_PC`.
3. Changed the fastENLOC comparison `slope[i]==-slope[i]` to the intended assignment `slope[i] <- -slope[i]`.
4. Fixed legacy ECLIPSER numbering, an extra gsMap subdirectory, broken HDL/PLACO/Figure 5 input links, and the PCa `Normal -> Healthy` label mismatch.
5. The PLACO loop generates `PCA_PSA`, consistent with the directories read by steps 31, 32, and 34; step 32 stops when files are missing rather than silently skipping them.
6. Continuous-PSA `06_prepare_ecaviar` processes continuous PSA only and excludes the three main-analysis traits.
7. Table 1 no longer uses TwoSampleMR, MRPRESSO, the two MR RDS files, or MR supplementary tables; it is built directly from eCAVIAR/fastENLOC results.
8. Total gene counts in Table 1 use CLPP>=0.01 or RCP>=0.10. BPH/PCa supplementary tables dynamically select the other two traits and do not generate self-colocalization columns.
9. Coordinates are retained for fastENLOC-only records in Table 1, preventing their loss during coordinate-based filtering. Display ranking uses the maximum single-method value after threshold normalization without changing inclusion thresholds.
10. Threshold sensitivity is documented as CLPP relaxed to 0.05 and RCP tightened to 0.50, with both thresholds shown in the output directory.
11. `Any_Tissue` in the scoring output was replaced by the accurate term `Non_Prostate`; colocalization and scDRS inputs are linked directly to final outputs.
12. Timestamps were removed from continuous-PSA LDSC/HDL output directories to ensure reproducible paths.
13. Personal absolute paths and dependencies on the active RStudio document were removed; R scripts end with `sessionInfo()`.
14. The COJO fallback EAF input is standardized to `main_analysis/rawdata/sumstats_rehab/result/`; the gsMap spatial transcriptomics input is standardized to `main_analysis/27_prepare_gsmap/ST/`.

## Static checks

The 67 pre-existing R scripts passed syntax-tree parsing, and the added subcluster merge-back R script passed quotation-mark and delimiter checks. All 9 shell scripts passed `bash -n`, and all 5 deduplicated embedded Python blocks passed syntax parsing. No broken R `source()` paths were found. Scans found no personal paths, unrelated inputs, MR dependencies, legacy directory names, non-English filenames, Chinese text, or key directory breaks. Continuous-PSA direction handling matches the confirmed orientation.

Static checks do not constitute numerical reproduction. The current environment lacks the complete original GWAS, GTEx, LD panels, and single-cell objects, and no claim is made that the full workflow has been executed successfully.

## Required local reruns or confirmations

1. **LAVA:** the uploaded version ran only PC-BPH, whereas the manuscript requires two pairs. Rerun with `BPH37.txt`, `PC37.txt`, and `PSA37.txt` and reconcile the outputs with the manuscript.
2. **fastENLOC:** the assignment error has been corrected. Rerun the affected preprocessing and downstream steps; numerical equivalence cannot be assumed in advance.
3. **Continuous-PSA direction:** swap `effect_allele` and `other_allele` according to the MVP file, with `EAF = 1 - EAF` and `beta = -beta`; rerun all continuous-PSA downstream analyses.
4. **Thirteen genes in the main score:** recalculate from the final inputs. If the number differs from 13, check the input versions first and do not add or remove genes manually.
5. **PLACO and Table 1:** after correcting directories, counts, and supplementary-table columns, rerun the affected tables and figures. Previous outputs must not be presented as outputs of the revised code.
