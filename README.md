# Cross-Trait Genetic Analysis of Prostate Traits

This repository accompanies the manuscript *Integrating genetic regulation and single-cell expression with GWAS prioritizes effector genes and cell types underlying PSA elevation in benign prostatic hyperplasia and prostate cancer*.

## 1. Repository scope

The analysis code is publicly available in this repository. The scripts have been deduplicated, renamed in English, sequentially numbered, reordered by dependencies, checked for input-output alignment, and statically validated.

Numerical reproduction requires the original GWAS, QTL, linkage disequilibrium reference, single-cell, and other input resources listed in `INPUT_MANIFEST.tsv`. Implemented corrections and analyses requiring local regeneration are documented in `CODE_AUDIT.md`.

The main analysis, continuous-PSA sensitivity analysis, and threshold sensitivity analysis write to separate directories and do not overwrite one another. Results must never be edited manually to match an earlier table.

## 2. Input GWAS datasets

| Analysis | Trait | Dataset | Sample size |
| --- | --- | --- | --- |
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

GWAS and QTL datasets were harmonized using GRCh38 chromosome, position, and allele identifiers. Effect/other allele orientation was harmonized according to the input requirements of each analysis. GWAS Z scores were calculated as `beta / se`, and QTL Z scores were calculated as `slope / slope_se`. The genome-build requirements and reference resources for individual analyses are specified in `PARAMETERS.md` and the corresponding scripts.

The 1000 Genomes European reference panel was used for eCAVIAR, fastENLOC, and conditional analyses. Additional analysis-specific linkage disequilibrium references are listed in `INPUT_MANIFEST.tsv`.

Colocalization analyses used a ±50-kb window around each lead SNP. The eCAVIAR settings were `causal = 2`, `rho_prob = 0.95`, and `f = 1`. Default priors were used for fastENLOC.

The primary gene-prioritization thresholds were CLPP ≥ 0.01 and RCP ≥ 0.10. The threshold sensitivity analysis used CLPP ≥ 0.05 and RCP ≥ 0.50, representing tightening of both thresholds.

## 5. Run order

Set the R working directory to the directory containing the script. Complete each numbered step before proceeding. For web-based steps, download the actual results and place them in the specified directories; do not substitute fictitious web automation.

### 5.1 Main analysis

1. `01_preprocess_gwas -> 02_harmonize_gwas`: format and quality-control the original GWAS files.
2. `03_qtl_enrich_eqtl / 04_qtl_enrich_sqtl -> 05_plot_qtl_enrich_eqtl / 06_plot_qtl_enrich_sqtl`: QTL enrichment and QQ plots.
3. `07_fuma_gwas -> 08_prepare_ecaviar -> 09_run_ecaviar -> 10_process_ecaviar`: FUMA locus processing and eCAVIAR.
4. `11_prepare_fastenloc -> 12_run_fastenloc -> 13_process_fastenloc`: fastENLOC.
5. `14_gene_enrich -> 15_process_gene_enrich -> 16_summarize_coloc`: GeneEnrich and colocalization summaries.
6. `17_prepare_coloc_plot` through `28_run_gsmap`: colocalization plots, xQTLbiolinks, HyPrColoc, ECLIPSER, COJO, and gsMap, as required by the manuscript analyses.
7. `29_ldsc` and `30_hdl`: genetic correlation analyses, which can be run independently after step 02.
8. `31_placo -> 32_prepare_placo_fuma`: PLACO and FUMA submission files.
9. `33_lava/01_format_lava_input.R -> 33_lava/02_run_lava_pairs.R`: local genetic correlations for PSA-BPH and PSA-PCa.
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

- `01_score_binary_psa.R`: main score with CLPP ≥ 0.01, RCP ≥ 0.10, and scDRS top 1000.
- `02_score_continuous_psa.R`: continuous-PSA sensitivity analysis using the same scoring rules.
- `03_score_threshold_sensitivity.R`: threshold sensitivity analysis with CLPP ≥ 0.05 and RCP ≥ 0.50.

The score directly reads the final outputs of main-analysis steps 10, 13, and 35 and the scDRS results. Continuous-PSA scoring reads sensitivity-analysis steps 08 and 11. Place the PLACO-MAGMA web outputs in `scoring/` as specified in the input manifest.

Within each trait pair, candidate genes are selected by ranking rather than by a statistical-significance cutoff. Among genes with a tissue/cell evidence score (`T_score`) of at least 2, the genes with the maximum total score are selected for the main-text candidate set. The union of the selected genes from the two trait pairs contains 13 genes in the current main analysis. The `Table2_Main_Text.csv` display rule of total score ≥7 is equivalent to this maximum-score selection in the current main result; 7 is not interpreted as a statistical-significance or causality threshold.

The complete score distributions are retained so that score and rank changes can be examined under the strict-colocalization and continuous-PSA analyses without redefining the main candidate set post hoc.

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

The total score therefore ranges from 0 to 8. The `cross` label and associated output columns refer specifically to non-prostate tissues. Prostate-tissue evidence is scored separately.

The scoring pool summarizes genes represented in the colocalization, Open4Gene, PLACO-MAGMA, and scDRS inputs. PLACO-MAGMA is used to describe pleiotropic genes and construct the scDRS input; it does not contribute points to the eight-component score. FUMA is used for locus definition and is not scored.

“Shared candidate effector genes” is an umbrella term for genes included in the multilayer evidence evaluation of the BPH-PSA or PCa-PSA trait pair. It denotes eligibility for statistical prioritization, not an experimentally established effector gene.

“Main-text prioritized genes” refers to the genes selected by the within-pair maximum-total-score rule described above. Where the manuscript uses “high-confidence genes” or “statistical-priority genes,” these labels refer only to this main-text display set. They do not denote a statistical confidence probability, independent clinical replication, or experimental proof of causality.

The total score is used for descriptive ranking. It is not an effect size, a statistical-significance measure, or a causal probability. eCAVIAR and fastENLOC use the same GWAS and QTL inputs, prostate-tissue evidence is a subset of the cross-tissue evidence, and scDRS is related to the upstream PLACO-MAGMA input. The eight components therefore must not be interpreted as eight independent statistical validations.

## 7. Interpretation boundaries

- Single-cell expression anchoring is descriptive. Large numbers of cells from the same patient are not treated as independent biological replicates for P-value calculation.
- The Hu healthy adult prostate atlas retains its 16 fine-grained cell-type annotations. It is used only to assess whether expression distributions are similar in the corresponding cellular compartments; the highest-expressing cell type does not redefine the original localization.
- Continuous-PSA analysis compares score and rank stability. It is not disease-stratified analysis or independent clinical validation.
- Colocalization, priority scores, and expression anchoring do not establish experimentally confirmed causal genes.

## 8. Version reporting

Analysis scripts record `sessionInfo()` directly or through shared entry points. Actual versions of command-line software must be recorded during the final runs. Uniform version numbers that cannot be verified from the uploaded scripts are not invented.
