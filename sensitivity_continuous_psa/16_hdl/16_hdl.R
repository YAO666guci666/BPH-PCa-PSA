library(QTLMR)
library(devtools)
library(data.table)

out_root <- "./HDL_PSA6319"

HDL_out <- file.path(out_root, "HDL")
HDL_L1_out <- file.path(out_root, "HDL_L", "1")
HDL_L2_out <- file.path(out_root, "HDL_L", "2")

dir.create(
  HDL_out,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  HDL_L1_out,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  HDL_L2_out,
  recursive = TRUE,
  showWarnings = FALSE
)

sink(
  file.path(out_root, "run_log.txt"),
  split = TRUE
)

cat("Working directory:", getwd(), "\n")
cat("Info Output directory:", normalizePath(out_root), "\n\n")

input_files <- c(
  "../../main_analysis/01_preprocess_gwas/BPH.txt",
  "../../main_analysis/01_preprocess_gwas/PCA.txt",
  "../01_preprocess_gwas/PSA_continuous.txt"
)

for (input_file in input_files) {
  if (!file.exists(input_file)) {
    stop(
      "Input file not found:",
      input_file,
      "\n Working directory:",
      getwd()
    )
  }

  cat("\n input file:", input_file, "\n")
  cat("file Info:", file.info(input_file)$size, "bytes\n")

  input_head <- fread(
    input_file,
    nrows = 3,
    data.table = FALSE
  )

  cat("Actual columns: \n")
  print(names(input_head))
  cat("Info 3 rows: \n")
  print(input_head)
}

dat_h2 <- HDL_h2(
  GWASfile = "../../main_analysis/01_preprocess_gwas/BPH.txt",
  LD.path = "./UKB_imputed_SVD_eigen99_extraction",
  Nref = 335265,
  eigen.cut = "automatic",
  intercept.output = FALSE,
  fill.missing.N = NULL,
  lim = exp(-18),
  save_name = "BPH",
  save_path = HDL_out
)

dat_h2 <- HDL_h2(
  GWASfile = "../../main_analysis/01_preprocess_gwas/PCA.txt",
  LD.path = "./UKB_imputed_SVD_eigen99_extraction",
  Nref = 335265,
  eigen.cut = "automatic",
  intercept.output = FALSE,
  fill.missing.N = NULL,
  lim = exp(-18),
  save_name = "PCA",
  save_path = HDL_out
)

dat_h2 <- HDL_h2(
  GWASfile = "../01_preprocess_gwas/PSA_continuous.txt",
  LD.path = "./UKB_imputed_SVD_eigen99_extraction",
  Nref = 335265,
  eigen.cut = "automatic",
  intercept.output = FALSE,
  fill.missing.N = NULL,
  lim = exp(-18),
  save_name = "PSA_6319",
  save_path = HDL_out
)

dat_rg <- HDL_rg_parallel(
  GWASfile = c(
    "../../main_analysis/01_preprocess_gwas/BPH.txt",
    "../01_preprocess_gwas/PSA_continuous.txt"
  ),
  LD.path = "./UKB_imputed_SVD_eigen99_extraction",
  Nref = 335265,
  N0 = NULL,
  numCores = 10,
  eigen.cut = "automatic",
  jackknife.df = FALSE,
  intercept.output = FALSE,
  fill.missing.N = NULL,
  lim = exp(-18),
  verbose = FALSE,
  save_name = "BPH_PSA_6319",
  save_path = HDL_out
)

dat_rg <- HDL_rg_parallel(
  GWASfile = c(
    "../../main_analysis/01_preprocess_gwas/PCA.txt",
    "../01_preprocess_gwas/PSA_continuous.txt"
  ),
  LD.path = "./UKB_imputed_SVD_eigen99_extraction",
  Nref = 335265,
  N0 = NULL,
  numCores = 10,
  eigen.cut = "automatic",
  jackknife.df = FALSE,
  intercept.output = FALSE,
  fill.missing.N = NULL,
  lim = exp(-18),
  verbose = FALSE,
  save_name = "PCA_PSA_6319",
  save_path = HDL_out
)

dat_L1 <- HDL_L(
  GWASfile = c(
    "../../main_analysis/01_preprocess_gwas/BPH.txt",
    "../01_preprocess_gwas/PSA_continuous.txt"
  ),
  Trait1name = "BPH",
  Trait2name = "PSA_6319",
  LD.path = "./UKB_REF/LD.path",
  bim.path = "./UKB_REF/bimfile",
  Nref = 335272,
  N0 = NULL,
  start_loc = NULL,
  end_loc = NULL,
  eigen.cut = 0.99,
  intercept.output = FALSE,
  fill.missing.N = NULL,
  lim = exp(-18),
  alpha = 0.05,
  save_name = "BPH_PSA_6319",
  save_path = HDL_L1_out,
  cores = 1
)

dat_L1 <- HDL_L(
  GWASfile = c(
    "../../main_analysis/01_preprocess_gwas/PCA.txt",
    "../01_preprocess_gwas/PSA_continuous.txt"
  ),
  Trait1name = "PCA",
  Trait2name = "PSA_6319",
  LD.path = "./UKB_REF/LD.path",
  bim.path = "./UKB_REF/bimfile",
  Nref = 335272,
  N0 = NULL,
  start_loc = NULL,
  end_loc = NULL,
  eigen.cut = 0.99,
  intercept.output = FALSE,
  fill.missing.N = NULL,
  lim = exp(-18),
  alpha = 0.05,
  save_name = "PCA_PSA_6319",
  save_path = HDL_L1_out,
  cores = 1
)

dat_L2 <- HDL_L(
  GWASfile = c(
    "../../main_analysis/01_preprocess_gwas/BPH.txt",
    "../01_preprocess_gwas/PSA_continuous.txt"
  ),
  Trait1name = "BPH",
  Trait2name = "PSA_6319",
  LD.path = "./UKB_PPP/LD",
  bim.path = "./UKB_PPP/bim",
  Nref = 335272,
  N0 = NULL,
  start_loc = NULL,
  end_loc = NULL,
  eigen.cut = 0.99,
  intercept.output = FALSE,
  fill.missing.N = NULL,
  lim = exp(-18),
  alpha = 0.05,
  save_name = "BPH_PSA_6319",
  save_path = HDL_L2_out,
  cores = 1
)

dat_L2 <- HDL_L(
  GWASfile = c(
    "../../main_analysis/01_preprocess_gwas/PCA.txt",
    "../01_preprocess_gwas/PSA_continuous.txt"
  ),
  Trait1name = "PCA",
  Trait2name = "PSA_6319",
  LD.path = "./UKB_PPP/LD",
  bim.path = "./UKB_PPP/bim",
  Nref = 335272,
  N0 = NULL,
  start_loc = NULL,
  end_loc = NULL,
  eigen.cut = 0.99,
  intercept.output = FALSE,
  fill.missing.N = NULL,
  lim = exp(-18),
  alpha = 0.05,
  save_name = "PCA_PSA_6319",
  save_path = HDL_L2_out,
  cores = 1
)

cat("\n Info rows Info. \n")
cat("Info Output directory:", normalizePath(out_root), "\n")
cat("Info Generating file: \n")
print(
  list.files(
    out_root,
    recursive = TRUE,
    full.names = TRUE
  )
)

sink()

sessionInfo()
