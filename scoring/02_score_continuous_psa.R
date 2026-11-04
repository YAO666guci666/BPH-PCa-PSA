rm(list = ls())
OUT_DIR <- "./scoring_continuous_psa_clpp0.01_rcp0.10"
TOP_N <- 1000
CLPP_CUT <- 0.01
RCP_CUT <- 0.10
EF_FILES <- list(
  BPH = c(cav = "../main_analysis/12_process_ecaviar/stable_eCAVIAR_finngen_BPH.sig.txt", fen = "../main_analysis/15_process_fastenloc/stable_finngen_BPH.fastenloc.sig.txt"),
  PCa = c(cav = "../main_analysis/12_process_ecaviar/stable_eCAVIAR_ieu-b-4809_PC.sig.txt", fen = "../main_analysis/15_process_fastenloc/stable_ieu-b-4809_PC.fastenloc.sig.txt"),
  PSA = c(cav = "../sensitivity_continuous_psa/10_process_ecaviar/stable_eCAVIAR_GCST90476319_PSA.sig.txt", fen = "../sensitivity_continuous_psa/13_process_fastenloc/stable_GCST90476319_PSA.fastenloc.sig.txt")
)
SCDRS_FILES <- c(BPH_PSA_6319 = "../scrna/scdrs/bph/scDRS_out/BPH_PSA_6319.bph_top1000_genes.txt", PCa_PSA_6319 = "../scrna/scdrs/pca/scDRS_out/PCA_PSA_6319.tumor_top1000_genes.txt")
PAIRS <- list(list(pair = "BPH_PSA_6319", T1 = "BPH", T2 = "PSA"), list(pair = "PCa_PSA_6319", T1 = "PCa", T2 = "PSA"))
O4G_FILE <- "../main_analysis/37_open4gene/Open4Gene.obj.res.txt"
EXPECTED_MAIN_GENES <- NA_integer_
CAPTION_TEXT <- "Table S. Candidate gene prioritization using continuous PSA (GCST90476319; Total Score >= 7).
Values in eCAVIAR and fastENLOC columns are presented as BPH / PSA / PCa and calculated from non-prostate tissues. A dash indicates that the corresponding trait did not reach the threshold or was not detected. Open4Gene indicates whether the gene was identified by Open4Gene analysis. scDRS correlation values are displayed only for genes ranked within the top 1000 of the corresponding single-cell dataset.

Scoring criteria are identical to the main analysis: eCAVIAR CLPP >= 0.01, fastENLOC RCP >= 0.10, Open4Gene, prostate colocalization and scDRS top 1000."
source("score_core.R")
