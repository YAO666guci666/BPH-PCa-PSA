rm(list = ls())

library(GenomicSEM)
library(data.table)
library(dplyr)
library(stringr)

out_dir <- "./LDSC_GCST90476319"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
sink(log_file, split = TRUE)

cat("Working directory:", getwd(), "\n")
cat("Info Output directory:", normalizePath(out_dir), "\n\n")

munge_data <- function(file, hm3, N, trait.name, out_dir, overwrite = FALSE) {
  ref <- fread(hm3, header = TRUE, data.table = FALSE)

  if (!"MAF" %in% colnames(file)) {
    if ("eaf" %in% colnames(file)) {
      file$MAF <- ifelse(file$eaf <= 0.5, file$eaf, 1 - file$eaf)
    }
  }

  if ("MAF" %in% colnames(file)) {
    file$MAF <- ifelse(file$MAF <= 0.5, file$MAF, 1 - file$MAF)
  }

  file$A1 <- factor(
    toupper(file$effect_allele),
    levels = c("A", "C", "G", "T")
  )
  file$A2 <- factor(
    toupper(file$other_allele),
    levels = c("A", "C", "G", "T")
  )

  file <- merge(ref, file, by = "SNP", all.x = FALSE, all.y = FALSE)

  file <- subset(file, !is.na(file$pval))
  file <- subset(file, !is.na(file$beta))

  file$beta <- ifelse(
    file$A1.x != file$A1.y & file$A1.x == file$A2.y,
    file$beta * -1,
    file$beta
  )

  file <- subset(
    file,
    !(file$A1.x != file$A1.y & file$A1.x != file$A2.y)
  )
  file <- subset(
    file,
    !(file$A2.x != file$A2.y & file$A2.x != file$A1.y)
  )

  file$Z <- sign(file$beta) * sqrt(
    qchisq(file$pval, 1, lower.tail = FALSE)
  )

  if ("N" %in% colnames(file)) {
    output <- cbind.data.frame(
      file$SNP,
      file$N,
      file$Z,
      file$A1.x,
      file$A2.x
    )
  } else {
    output <- cbind.data.frame(
      file$SNP,
      N,
      file$Z,
      file$A1.x,
      file$A2.x
    )
  }

  colnames(output) <- c("SNP", "N", "Z", "A1", "A2")

  trait.name <- str_replace_all(trait.name, fixed(" "), "")

  sumstats_file <- file.path(
    out_dir,
    paste0(trait.name, ".sumstats")
  )

  fwrite(
    x = output,
    file = sumstats_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  R.utils::gzip(
    sumstats_file,
    overwrite = overwrite
  )

  cat(
    trait.name,
    "After preprocessing SNP count:",
    nrow(output),
    "\n"
  )

  return()
}

traits_info <- list(
  BPH = list(
    file = "../../main_analysis/02_harmonize_gwas/finngen_BPH_hg38_qc.txt",
    N = 194710,
    name = "BPH",
    sample.prev = 0.2113,
    population.prev = 0.40
  ),
  PC = list(
    file = "../../main_analysis/02_harmonize_gwas/ieu-b-4809_PC_hg38_qc.txt",
    N = 182625,
    name = "PC",
    sample.prev = 0.0500,
    population.prev = 0.08
  ),
  PSA_cont = list(
    file = "../02_harmonize_gwas/GCST90476319_PSA_hg38_qc.txt",
    N = 254502,
    name = "PSA_cont",
    sample.prev = NA,
    population.prev = NA
  )
)

hm3_path <- "eur_w_ld_chr/w_hm3.snplist"

required_columns <- c(
  "SNP",
  "effect_allele",
  "other_allele",
  "beta",
  "pval",
  "eaf"
)

for (trait_name in names(traits_info)) {
  trait <- traits_info[[trait_name]]

  cat("\nProcessing", trait$name, "...\n")

  if (!file.exists(trait$file)) {
    stop(
      "Input file not found:",
      trait$file,
      "\n Working directory:",
      getwd()
    )
  }

  trait_file <- fread(
    trait$file,
    data.table = FALSE
  )

  cat("file:", trait$file, "\n")
  cat("rows, columns Info:", nrow(trait_file), ncol(trait_file), "\n")
  cat("Actual columns: \n")
  print(names(trait_file))

  missing_columns <- setdiff(
    required_columns,
    names(trait_file)
  )

  if (length(missing_columns) > 0) {
    stop(
      "file",
      trait$file,
      "Missing the following columns:",
      paste(missing_columns, collapse = ", "),
      "\n Check the actual columns printed above, do not continue."
    )
  }

  cat("Info 3 rows Info columns: \n")
  print(
    head(
      trait_file[, required_columns],
      3
    )
  )

  munge_data(
    trait_file,
    hm3 = hm3_path,
    N = trait$N,
    trait.name = trait$name,
    out_dir = out_dir
  )
}

trait_pairs <- list(
  c("BPH", "PSA_cont"),
  c("PC", "PSA_cont")
)

all_results <- list()

for (pair in trait_pairs) {
  trait1_name <- pair[1]
  trait2_name <- pair[2]

  cat(
    "\nRunning LDSC for",
    trait1_name,
    "and",
    trait2_name,
    "...\n"
  )

  traits <- c(
    file.path(
      out_dir,
      paste0(traits_info[[trait1_name]]$name, ".sumstats.gz")
    ),
    file.path(
      out_dir,
      paste0(traits_info[[trait2_name]]$name, ".sumstats.gz")
    )
  )

  sample.prev <- c(
    traits_info[[trait1_name]]$sample.prev,
    traits_info[[trait2_name]]$sample.prev
  )

  population.prev <- c(
    traits_info[[trait1_name]]$population.prev,
    traits_info[[trait2_name]]$population.prev
  )

  trait.names <- c(
    traits_info[[trait1_name]]$name,
    traits_info[[trait2_name]]$name
  )

  tryCatch({
    LDSCoutput <- GenomicSEM::ldsc(
      traits = traits,
      sample.prev = sample.prev,
      population.prev = population.prev,
      ld = "eur_w_ld_chr/",
      wld = "eur_w_ld_chr/",
      trait.names = trait.names,
      stand = TRUE
    )

    h2_exp <- as.numeric(LDSCoutput$S[1, 1])
    h2_exp_SE <- sqrt(LDSCoutput$V[1, 1])
    gcov_int <- as.numeric(LDSCoutput$I[1, 2])
    int_exp <- LDSCoutput$I[1, 1]
    int_out <- LDSCoutput$I[2, 2]
    h2_out <- as.numeric(LDSCoutput$S[2, 2])
    h2_out_SE <- sqrt(LDSCoutput$V[3, 3])
    rgcov <- as.numeric(LDSCoutput$S[1, 2])

    rg <- as.numeric(
      LDSCoutput$S[1, 2] /
        sqrt(
          LDSCoutput$S[1, 1] *
            LDSCoutput$S[2, 2]
        )
    )

    k <- nrow(LDSCoutput$S)
    SE <- matrix(0, k, k)
    SE[lower.tri(SE, diag = TRUE)] <- sqrt(
      diag(LDSCoutput$V_Stand)
    )

    rg_SE <- SE[2, 1]
    rg_P <- 2 * pnorm(
      q = abs(rg / rg_SE),
      lower.tail = FALSE
    )

    rgcov_SE <- sqrt(LDSCoutput$V[2, 2])

    ldsc_result <- data.frame(
      trait1 = trait1_name,
      trait2 = trait2_name,
      h2_trait1 = h2_exp,
      h2_trait1_se = h2_exp_SE,
      h2_trait2 = h2_out,
      h2_trait2_se = h2_out_SE,
      rg = rg,
      rg_se = rg_SE,
      rg_p = rg_P,
      rgcov = rgcov,
      rgcov_se = rgcov_SE,
      gcov_int = gcov_int,
      int_trait1 = int_exp,
      int_trait2 = int_out
    )

    all_results[[paste(
      trait1_name,
      trait2_name,
      sep = "_"
    )]] <- ldsc_result

    output_file <- file.path(
      out_dir,
      paste0(
        "LDSC_",
        trait1_name,
        "_",
        trait2_name,
        "_result.txt"
      )
    )

    write.table(
      ldsc_result,
      output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    cat("Results saved to", output_file, "\n")

    cat("Info results: \n")
    print(
      ldsc_result[, c(
        "trait1",
        "trait2",
        "rg",
        "rg_se",
        "rg_p",
        "gcov_int"
      )]
    )

  }, error = function(e) {
    cat(
      "Error in LDSC for",
      trait1_name,
      "and",
      trait2_name,
      ":",
      e$message,
      "\n"
    )
  })
}

if (length(all_results) > 0) {
  combined_results <- do.call(
    rbind,
    all_results
  )

  combined_file <- file.path(
    out_dir,
    "all_LDSC_results_GCST90476319.txt"
  )

  write.table(
    combined_results,
    combined_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  mydata <- combined_results %>%
    dplyr::select(
      trait1,
      trait2,
      rg,
      rg_se,
      rg_p
    )

  mydata$rg_lower <- mydata$rg - 1.96 * mydata$rg_se
  mydata$rg_upper <- mydata$rg + 1.96 * mydata$rg_se

  mydata$rg <- round(mydata$rg, 3)
  mydata$rg_se <- round(mydata$rg_se, 3)

  mydata$rg_p <- formatC(
    mydata$rg_p,
    format = "e",
    digits = 2
  )

  mydata$pair <- paste(
    mydata$trait1,
    "vs",
    mydata$trait2
  )

  mydata <- mydata %>%
    dplyr::select(
      pair,
      rg,
      rg_se,
      rg_p,
      rg_lower,
      rg_upper
    )

  mydata$` ` <- paste(
    rep(" ", 20),
    collapse = " "
  )

  library(forestploter)

  forest_file <- file.path(
    out_dir,
    "LDSC_forest_plot_GCST90476319.pdf"
  )

  pdf(
    forest_file,
    width = 10,
    height = 8
  )

  forest(
    mydata[, c(1, 2, 7, 5, 6)],
    est = mydata$rg,
    lower = mydata$rg_lower,
    upper = mydata$rg_upper,
    sizes = 0.5,
    ci_column = 3,
    ref_line = 0,
    xlim = c(-1, 1),
    arrow_lab = c(
      "Negative correlation",
      "Positive correlation"
    ),
    xlab = "Genetic Correlation (rg)"
  )

  dev.off()

  cat("\nAll analyses completed!\n")
  cat("Info results Directory:", normalizePath(out_dir), "\n")
  cat("Merging results:", combined_file, "\n")
  cat("Info:", forest_file, "\n")

} else {
  cat("\nNo successful LDSC analyses completed.\n")
}

sink()

sessionInfo()
