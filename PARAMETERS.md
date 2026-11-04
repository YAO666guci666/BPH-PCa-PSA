# Parameters and Reference Resources

## Data harmonization and quality control

- Main downstream analyses use GRCh38; liftOver is not repeated when the original input is already GRCh38.
- rsID, chromosome position, and effect/other alleles are used for matching. Any direction reversal must update beta and EAF together.
- MHC region: chr6:25,000,000-35,000,000.
- MAF filter: retain 0.01 < EAF < 0.99.
- For continuous PSA, `effect_allele` and `other_allele` are swapped according to the MVP file, with `EAF = 1 - EAF` and `beta = -beta` updated simultaneously. The sensitivity analysis must be rerun using this orientation.

## LD references

| Analysis | Reference |
|---|---|
| eCAVIAR, fastENLOC, COJO, locusZoom | 1000 Genomes EUR, `g1000_eur` |
| HyPrColoc | `./1kg.v3/EUR` |
| LDSC | European HM3 variants and `eur_w_ld_chr` |
| Genome-wide HDL | UKB reference, `Nref=335265` |
| Local HDL | UKB reference, `Nref=335272`, `eigen.cut=0.99` |
| LAVA | `EUR` and the predefined locus file |

## Colocalization and loci

- Each colocalization window extends 50 kb upstream and downstream of the lead SNP.
- GWAS and QTL data are harmonized by variant and allele. At least five variants are required in the LD matrix.
- eCAVIAR: `causal=2`, `rho_prob=0.95`, `f=1`; main-score threshold CLPP>=0.01.
- fastENLOC: default colocalization priors, `thread=8`; main-score threshold RCP>=0.10.
- HyPrColoc: `KLK3/rs17632542`, `JAZF1/rs10486567`, `CLPTM1L/rs414965`, and `NKX3-1/rs995432`; `binary.outcomes=c(TRUE,TRUE,TRUE)`, `ind.traits=TRUE`; LD r2=0.1.
- PLACO: `p_threshold=5e-8`, `filter_Zscores=TRUE`, `FDR_method="fdr"`, `decor=FALSE`, `cores=4`.
- FUMA sample sizes after PLACO: BPH-PSA=596,163, BPH-PCa=377,335, and PCa-PSA=584,078; 1000G EUR, GRCh38.

## Enrichment analyses

- GTEx v8 eQTL/sQTL data across 49 tissues.
- GeneEnrich: `min_permutations=1000`, `max_permutations=100000`, `is_HLA_remove=TRUE`.
- GENCODE v26.
- Empirical P-values are retained. BH-FDR is calculated within each trait-by-QTL-type-by-pathway-database stratum; empirical P-values do not replace cross-pathway BH-FDR.

## Genetic correlation and LAVA

- LDSC: BPH `N=194710, sample.prev=0.2113, population.prev=0.40`; PSA `N=401453, sample.prev=0.1355, population.prev=0.15`; PCa `N=182625, sample.prev=0.0500, population.prev=0.08`; `stand=T`.
- Continuous PSA: `N=254502`, with both prevalence values set to NA.
- LAVA: univariate threshold `0.05/nrow(loci)`; bivariate results are summarized using BH-FDR.

## Single-cell analysis and scDRS

- Seurat, Harmony, and scDblFinder. Dataset-specific clustering resolutions and manual annotation mappings remain in their original scripts and are not forced to be identical.
- scDblFinder: `dbr=0.08`, `nfeatures=1000`.
- scDRS: gene sets are generated from MAGMA Z scores with `n-max=1000`; scoring uses `--n-ctrl 1000 --weight-opt vs`.
- Cell-level significance summaries use `mc_pval<0.1`. Cell-type group associations are reported separately using `assoc_mcp` and its BH-FDR; these two quantities are not interchanged.
- Top 1000: within disease-group cells, genome-wide expression is correlated with normalized scDRS scores, and the top 1000 genes are selected in descending order of correlation coefficient.

## Eight-component score

Each component contributes one point, giving a total score from 0 to 8:

1. Maximum non-prostate eCAVIAR CLPP>=0.01 for BPH.
2. Maximum non-prostate eCAVIAR CLPP>=0.01 for PSA.
3. Maximum non-prostate fastENLOC RCP>=0.10 for BPH.
4. Maximum non-prostate fastENLOC RCP>=0.10 for PSA.
5. Open4Gene hit.
6. Prostate-tissue CLPP>=0.01 or RCP>=0.10 for BPH.
7. Prostate-tissue CLPP>=0.01 or RCP>=0.10 for PSA.
8. Membership in the scDRS top 1000 for the corresponding trait pair.

Threshold sensitivity analysis separately applies CLPP>=0.05 and RCP>=0.50. This mixed scheme relaxes one threshold and tightens the other and does not overwrite the main score.

Genes with Total>=7 in either trait pair enter the main table. This rule is applied uniformly by code; genes must not be added or removed manually. `Candidate` separately denotes genes with at least two points in the T dimension and the joint highest Total within a trait pair; it is not the same field as the Total>=7 main-table inclusion rule.
