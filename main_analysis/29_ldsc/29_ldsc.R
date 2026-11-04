library(GenomicSEM)
library(data.table)
library(dplyr)
library(stringr)

munge_data <- function(file, hm3, N, trait.name, overwrite = TRUE) {
  ref <- fread(hm3, header = T, data.table = F)

  if (!"MAF" %in% colnames(file)) {
    if ("eaf" %in% colnames(file)) {
      file$MAF <- ifelse(file$eaf <= 0.5, file$eaf, 1 - file$eaf)
    }
  }

  if ("MAF" %in% colnames(file)) {
    file$MAF <- ifelse(file$MAF <= 0.5, file$MAF, (1 - file$MAF))
  }

  file$A1 <- factor(toupper(file$effect_allele), c("A", "C", "G", "T"))
  file$A2 <- factor(toupper(file$other_allele), c("A", "C", "G", "T"))

  file <- merge(ref, file, by = "SNP", all.x = F, all.y = F)

  file <- subset(file, !(is.na(file$pval)))
  file <- subset(file, !(is.na(file$beta)))

  file$beta <- ifelse(file$A1.x != (file$A1.y) & file$A1.x == (file$A2.y),
                      file$beta * -1, file$beta)

  file <- subset(file, !(file$A1.x != (file$A1.y) & file$A1.x != (file$A2.y)))
  file <- subset(file, !(file$A2.x != (file$A2.y) & file$A2.x != (file$A1.y)))

  file$Z <- sign(file$beta) * sqrt(qchisq(file$pval, 1, lower.tail = F))

  if ("N" %in% colnames(file)) {
    output <- cbind.data.frame(file$SNP, file$N, file$Z,
                               file$A1.x, file$A2.x)
  } else {
    output <- cbind.data.frame(file$SNP, N, file$Z, file$A1.x,
                               file$A2.x)
  }

  colnames(output) <- c("SNP", "N", "Z", "A1", "A2")

  trait.name <- str_replace_all(trait.name, fixed(" "), "")

  fwrite(x = output, file = paste0(trait.name, ".sumstats"),
         sep = "\t", quote = FALSE, row.names = F)
  R.utils::gzip(paste0(trait.name, ".sumstats"), overwrite = overwrite)

  return()
}

traits_info <- list(
  BPH = list(
    file = "../02_harmonize_gwas/finngen_BPH_hg38_qc.txt",
    N = 194710,
    name = "BPH",
    sample.prev = 0.2113,
    population.prev = 0.40
  ),
  PSA = list(
    file = "../02_harmonize_gwas/GCST90476277_PSA_hg38_qc.txt",
    N = 401453,
    name = "PSA",
    sample.prev = 0.1355,
    population.prev = 0.15
  ),
  PC = list(
    file = "../02_harmonize_gwas/ieu-b-4809_PC_hg38_qc.txt",
    N = 182625,
    name = "PC",
    sample.prev = 0.0500,
    population.prev = 0.08
  )
)

hm3_path <- "eur_w_ld_chr/w_hm3.snplist"

for (trait_name in names(traits_info)) {
  trait <- traits_info[[trait_name]]

  cat("Processing", trait$name, "...\n")
  trait_file <- fread(trait$file, data.table = F)

  trait_file <- trait_file %>%
    dplyr::rename(
      SNP = SNP,
      effect_allele = effect_allele,
      other_allele = other_allele,
      beta = beta,
      pval = pval,
      eaf = eaf
    )

  munge_data(trait_file, hm3 = hm3_path, N = trait$N, trait.name = trait$name)
}

trait_pairs <- combn(names(traits_info), 2, simplify = FALSE)

all_results <- list()

for (pair in trait_pairs) {
  trait1_name <- pair[1]
  trait2_name <- pair[2]

  cat("\nRunning LDSC for", trait1_name, "and", trait2_name, "...\n")

  traits <- c(
    paste0(traits_info[[trait1_name]]$name, ".sumstats.gz"),
    paste0(traits_info[[trait2_name]]$name, ".sumstats.gz")
  )

  sample.prev <- c(
    traits_info[[trait1_name]]$sample.prev,
    traits_info[[trait2_name]]$sample.prev
  )
  population.prev <- c(
    traits_info[[trait1_name]]$population.prev,
    traits_info[[trait2_name]]$population.prev
  )
  trait.names <- c(traits_info[[trait1_name]]$name, traits_info[[trait2_name]]$name)

  tryCatch({
    LDSCoutput <- GenomicSEM::ldsc(
      traits = traits,
      sample.prev = sample.prev,
      population.prev = population.prev,
      ld = "eur_w_ld_chr/",
      wld = "eur_w_ld_chr/",
      trait.names = trait.names,
      stand = T
    )

    h2_exp = as.numeric(LDSCoutput$S[1,1])
    h2_exp_SE = sqrt(LDSCoutput$V[1,1])
    gcov_int = as.numeric(LDSCoutput$I[1,2])
    int_exp = LDSCoutput$I[1,1]
    int_out = LDSCoutput$I[2,2]
    h2_out = as.numeric(LDSCoutput$S[2,2])
    h2_out_SE = sqrt(LDSCoutput$V[3,3])
    rgcov = as.numeric(LDSCoutput$S[1,2])
    rg = as.numeric(LDSCoutput$S[1,2]/sqrt(LDSCoutput$S[1,1]*LDSCoutput$S[2,2]))

    k <- nrow(LDSCoutput$S)
    SE <- matrix(0, k, k)
    SE[lower.tri(SE, diag = TRUE)] <- sqrt(diag(LDSCoutput$V_Stand))
    rg_SE <- SE[2,1]
    rg_P <- 2 * pnorm(q = abs(rg/rg_SE), lower.tail = FALSE)
    rgcov_SE <- sqrt(LDSCoutput$V[2,2])

    ldsc_result <- data.frame(
      "trait1" = trait1_name,
      "trait2" = trait2_name,
      "h2_trait1" = h2_exp,
      "h2_trait1_se" = h2_exp_SE,
      "h2_trait2" = h2_out,
      "h2_trait2_se" = h2_out_SE,
      "rg" = rg,
      "rg_se" = rg_SE,
      "rg_p" = rg_P,
      "rgcov" = rgcov,
      "rgcov_se" = rgcov_SE,
      "gcov_int" = gcov_int,
      "int_trait1" = int_exp,
      "int_trait2" = int_out
    )

    all_results[[paste(trait1_name, trait2_name, sep = "_")]] <- ldsc_result

    output_file <- paste0("LDSC_", trait1_name, "_", trait2_name, "_result.txt")
    write.table(ldsc_result, output_file, sep = "\t", quote = F, row.names = F, col.names = T)

    cat("Results saved to", output_file, "\n")

  }, error = function(e) {
    cat("Error in LDSC for", trait1_name, "and", trait2_name, ":", e$message, "\n")
  })
}

if (length(all_results) > 0) {
  combined_results <- do.call(rbind, all_results)

  write.table(combined_results, "all_LDSC_results.txt", sep = "\t", quote = F, row.names = F, col.names = T)

  mydata <- combined_results %>%
    dplyr::select(trait1, trait2, rg, rg_se, rg_p)

  mydata$rg_lower <- mydata$rg - 1.96 * mydata$rg_se
  mydata$rg_upper <- mydata$rg + 1.96 * mydata$rg_se

  mydata$rg <- round(mydata$rg, 3)
  mydata$rg_se <- round(mydata$rg_se, 3)
  mydata$rg_p <- round(mydata$rg_p, 3)

  mydata$pair <- paste(mydata$trait1, "vs", mydata$trait2)

  mydata <- mydata %>%
    dplyr::select(pair, rg, rg_se, rg_p, rg_lower, rg_upper)

  mydata$` ` <- paste(rep(" ", 20), collapse = " ")

  library(forestploter)

  pdf("LDSC_forest_plot_all_pairs.pdf", width = 10, height = 8)

  forest(mydata[, c(1, 2, 7, 5, 6)],
         est = mydata$rg,
         lower = mydata$rg_lower,
         upper = mydata$rg_upper,
         sizes = 0.5,
         ci_column = 3,
         ref_line = 0,
         xlim = c(-1, 1),
         arrow_lab = c("Negative correlation", "Positive correlation"),
         xlab = "Genetic Correlation (rg)"

  )

  dev.off()

  cat("\nAll analyses completed! Results saved to:\n")
  cat("- all_LDSC_results.txt (combined results)\n")
  cat("- LDSC_forest_plot_all_pairs.pdf (forest plot)\n")
  cat("- Individual result files for each trait pair\n")
} else {
  cat("\nNo successful LDSC analyses completed.\n")
}

sessionInfo()
