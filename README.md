# Cross-Trait Genetic Analysis of Prostate Traits

This repository accompanies the manuscript *Integrating genetic regulation and single-cell expression with GWAS prioritizes effector genes and cell types underlying PSA elevation in benign prostatic hyperplasia and prostate cancer*.

## 1. Repository scope

The analysis code is publicly available in this repository. The scripts have been deduplicated, renamed in English, sequentially numbered, reordered by dependencies, checked for input-output alignment, and statically validated.

Numerical reproduction requires the original GWAS, QTL, linkage disequilibrium reference, single-cell, and other input resources listed in `INPUT_MANIFEST.tsv`. Implemented corrections and analyses requiring local regeneration are documented in `CODE_AUDIT.md`.

The main analysis, continuous-PSA sensitivity analysis, and threshold sensitivity analysis write to separate directories and do not overwrite one another. Results must never be edited manually to match an earlier table.

## 2. Input GWAS datasets

| Analysis | Trait | Dataset | Sample size |
|---|---|---|---|
| Main | BPH | FinnGen R10, binary | 41,137 cases and 153,573 controls; N=194,710 |
| Main | PCa | IEU OpenGWAS `ieu-b-4809`, binary | 9,132 cases and 173,493 controls; N=182,625 |
| Main | Elevated PSA | GWAS Catalog `GCST90476277`, binary | 54,381 cases and 347,072 controls; N=401,453 |
| Sensitivity | Continuous PSA | GWAS Catalog `GCST90476319` | N=254,502 |

`BPH.txt`, `PSA.txt`, `PCA.txt`, and `PSA_continuous.txt` are downstream aliases of formatted original GWAS files, not MTAG joint statistics. This repository does not contain an MTAG analysis.

## 3. Directory structure and shared code

- `main_analysis/`: main analyses of the three binary traits, numbered continuously from `01_preprocess_gwas` to `36_table1_colocalization`.
- `sensitivity_continuous_psa/`: continuous-PSA preprocessing and shared-pipeline entry points, numbered from `01_preprocess_gwas` to `16_hdl`. Steps 02-14 call the corresponding main-analysis scripts with different trait inputs rather than duplicating the full code.
- `single_cell/`: BPH and PCa download, quality control, Harmony integration, broad annotation, epithelial subclustering, CopyKAT, subcluster merge-back, and expression anchoring.
- `scrna/scdrs/`: MAGMA conversion, scDRS scoring, cell-type statistics, and top-1000 export. Shared BPH/PCa matrix conversion, statistics, and top-1000 code are reused. Step 02 remains dataset-specific because the cell sets and downstream settings differ.
- `scoring/`: three short entry scripts for the main score, continuous-PSA score, and threshold sensitivity score. All calculations are centralized in `score_core.R`, which is not run directly.
- `PARAMETERS.md`: analysis parameters and reference panels.
- `INPUT_MANIFEST.tsv`: required raw data, reference resources, and web-analysis outputs.
- `CODE_AUDIT.md`: implemented corrections, static checks, and items requiring local regeneration.

## 4. Key analysis settings

GWAS and QTL datasets were harmonized using variant identifiers and effect/other alleles according to the input requirements of each analysis. The genome-build requirements and reference resources for individual analyses are specified in `PARAMETERS.md` and the corresponding scripts.

The 1000 Genomes European reference panel was used for eCAVIAR, fastENLOC, and conditional analyses. Additional analysis-specific linkage disequilibrium references are listed in `INPUT_MANIFEST.tsv`.

Colocalization analyses used a ±50-kb window around each lead SNP. The eCAVIAR settings were `causal = 2`, `rho_prob = 0.95`, and `f = 1`. Default priors were used for fastENLOC.

The primary gene-prioritization thresholds were CLPP ≥ 0.01 and RCP ≥ 0.10. The threshold sensitivity analysis used CLPP ≥ 0.05 and RCP ≥ 0.50, representing relaxation of the CLPP threshold and tightening of the RCP threshold rather than a uniform relaxation.

## 5. Run order

Set the R working directory to the directory containing the script. Complete each numbered step before proceeding. For web-based steps, download the actual results and place them in the specified directories; do not substitute fictitious web automation.

### 5.1 Main analysis

1. `01_preprocess_gwas -&gt; 02_harmonize_gwas`: format and quality-control the original GWAS files.
2. `03_qtl_enrich_eqtl / 04_qtl_enrich_sqtl -&gt; 05_plot_qtl_enrich_eqtl / 06_plot_qtl_enrich_sqtl`: QTL enrichment and QQ plots.
3. `07_fuma_gwas -&gt; 08_prepare_ecaviar -&gt; 09_run_ecaviar -&gt; 10_process_ecaviar`: FUMA locus processing and eCAVIAR.
4. `11_prepare_fastenloc -&gt; 12_run_fastenloc -&gt; 13_process_fastenloc`: fastENLOC.
5. `14_gene_enrich -&gt; 15_process_gene_enrich -&gt; 16_summarize_coloc`: GeneEnrich and colocalization summaries.
6. `17_prepare_coloc_plot` through `28_run_gsmap`: colocalization plots, xQTLbiolinks, HyPrColoc, ECLIPSER, COJO, and gsMap, as required by the manuscript analyses.
7. `29_ldsc` and `30_hdl`: genetic correlation analyses, which can be run independently after step 02.
8. `31_placo -&gt; 32_prepare_placo_fuma`: PLACO and FUMA submission files.
9. `33_lava/01_format_lava_input.R -&gt; 33_lava/02_run_lava_pairs.R`: local genetic correlations for PSA-BPH and PSA-PCa.
10. `34_placo_gene_figure`: integrate Figure 5 after the final outputs of steps 31, 32, and 33 are available.
11. `35_open4gene`: Open4Gene analysis.
12. `36_table1_colocalization`: generate Table 1 after steps 10 and 13; use colocalization results only and do not read MR objects.

### 5.2 Continuous-PSA sensitivity analysis

Run `01_preprocess_gwas`, followed by the shared-pipeline entry points in steps 02-14. Step `06_prepare_ecaviar` processes continuous PSA only and does not repeat the three main-analysis traits. Steps `15_ldsc` and `16_hdl` write to the fixed directories `LDSC_GCST90476319` and `HDL_PSA6319`, respectively, without timestamps.

### 5.3 Single-cell analysis and scDRS

1. Run the numbered scripts in `single_cell/bph` and `single_cell/pca`.
2. In each of these two directories, run `source("../06_merge_subclusters.R")` to generate `newBPH.rds` or `newPCa.rds` and `cell_type_sub_stats.csv`.
3. In `single_cell/`, run `07_expression_anchor_13genes.R` to read the merged objects from both datasets.
4. Run steps 01 and 02 in `scrna/scdrs/magma`.
5. Run steps 01-04 separately in `scrna/scdrs/bph` and `scrna/scdrs/pca`; each directory reads and writes only its own dataset files.

The main trait-pair names are `BPH_PSA` and `PCa_PSA`. The continuous-PSA names are `BPH_PSA_6319` and `PCA_PSA_6319`.

### 5.4 Scoring

Run one of the following entry scripts in `scoring/`:

- `01_score_binary_psa.R`: main score with CLPP ≥ 0.01, RCP ≥ 0.10, scDRS top 1000, and Total ≥ 7 in either trait pair.
- `02_score_continuous_psa.R`: continuous-PSA sensitivity analysis using the same scoring rules.
- `03_score_threshold_sensitivity.R`: mixed threshold sensitivity analysis with CLPP ≥ 0.05 and RCP ≥ 0.50.

The score directly reads the final outputs of main-analysis steps 10, 13, and 35 and the scDRS results. Continuous-PSA scoring reads sensitivity-analysis steps 08 and 11. Place the PLACO-MAGMA web outputs in `scoring/` as specified in the input manifest.

## 6. Gene prioritization and terminology

Eight evidence components are equally weighted, with one point assigned to each component:

1. eCAVIAR evidence for the first trait in non-prostate tissues.
2. eCAVIAR evidence for the second trait in non-prostate tissues.
3. fastENLOC evidence for the first trait in non-prostate tissues.
4. fastENLOC evidence for the second trait in non-prostate tissues.
5. Open4Gene evidence.
6. Prostate-tissue colocalization evidence for the first trait.
7. Prostate-tissue colocalization evidence for the second trait.
8. scDRS evidence for the corresponding trait pair.

The total score therefore ranges from 0 to 8. The `cross` label and associated output columns refer specifically to non-prostate tissues.

“Shared candidate effector genes” is used as an umbrella term for genes included in the multilayer evidence evaluation of the BPH-PSA or PCa-PSA trait pair to assess their potential involvement in shared genetic regulation.

“High-priority candidate genes” refer to the subset that meets the prespecified prioritization threshold. A gene is defined as high-priority when its total score is ≥7 in at least one of the two trait pairs.

These terms distinguish the broader candidate set from the statistically prioritized subset. The total score represents descriptive statistical prioritization, not an effect size or experimental evidence of causality.

## 7. Interpretation boundaries

- Single-cell expression anchoring is descriptive. Large numbers of cells from the same patient are not treated as independent biological replicates for P-value calculation.
- The Hu adult normal prostate atlas is used only to assess whether expression distributions are similar in the corresponding cellular compartments; the highest-expressing cell type does not redefine the original localization.
- Colocalization, priority scores, and expression anchoring do not establish experimentally confirmed causal genes.

## 8. Version reporting

Analysis scripts record `sessionInfo()` directly or through shared entry points. Actual versions of command-line software must be recorded during the final runs. Uniform version numbers that cannot be verified from the uploaded scripts are not invented.
