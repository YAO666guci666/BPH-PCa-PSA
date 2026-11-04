rm(list = ls())
OUT_DIR <- "./scoring_binary_psa_clpp0.01_rcp0.10"
TOP_N <- 1000
CLPP_CUT <- 0.01
RCP_CUT <- 0.10
EF_FILES <- list(
  BPH = c(cav = "../main_analysis/12_process_ecaviar/stable_eCAVIAR_finngen_BPH.sig.txt", fen = "../main_analysis/15_process_fastenloc/stable_finngen_BPH.fastenloc.sig.txt"),
  PCa = c(cav = "../main_analysis/12_process_ecaviar/stable_eCAVIAR_ieu-b-4809_PC.sig.txt", fen = "../main_analysis/15_process_fastenloc/stable_ieu-b-4809_PC.fastenloc.sig.txt"),
  PSA = c(cav = "../main_analysis/12_process_ecaviar/stable_eCAVIAR_GCST90476277_PSA.sig.txt", fen = "../main_analysis/15_process_fastenloc/stable_GCST90476277_PSA.fastenloc.sig.txt")
)
SCDRS_FILES <- c(BPH_PSA = "../scrna/scdrs/bph/scDRS_out/BPH_PSA.bph_top1000_genes.txt", PCa_PSA = "../scrna/scdrs/pca/scDRS_out/PCa_PSA.tumor_top1000_genes.txt")
PAIRS <- list(list(pair = "BPH_PSA", T1 = "BPH", T2 = "PSA"), list(pair = "PCa_PSA", T1 = "PCa", T2 = "PSA"))
O4G_FILE <- "../main_analysis/37_open4gene/Open4Gene.obj.res.txt"
EXPECTED_MAIN_GENES <- 13L
CAPTION_TEXT <- "Table 2. Candidate gene prioritization using binary PSA (Total Score >= 7).
Values in eCAVIAR and fastENLOC columns are presented as BPH / PSA / PCa and calculated from non-prostate tissues. A dash indicates that the corresponding trait did not reach the threshold or was not detected. Open4Gene indicates whether the gene was identified by Open4Gene analysis. scDRS correlation values are displayed only for genes ranked within the top 1000 of the corresponding single-cell dataset.

Scoring criteria: eCAVIAR CLPP >= 0.01 and fastENLOC RCP >= 0.10 in non-prostate tissues (+1 per trait and method); Open4Gene (+1); prostate colocalization using the same thresholds (+1 per trait); scDRS top 1000 (+1 per pair)."
source("score_core.R")
