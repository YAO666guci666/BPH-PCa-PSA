rm(list = ls())
sink("run_log.txt", split = TRUE)
on.exit(sink(), add = TRUE)
library(QTLMR)
library(data.table)
library(hyprcoloc)

PLINK <- "/home/train/miniconda3/bin/plink"
REF <- "./1kg.v3/EUR"

gwas_files <- c(
  BPH = "../02_harmonize_gwas/finngen_BPH_hg38_qc.txt",
  PCa = "../02_harmonize_gwas/ieu-b-4809_PC_hg38_qc.txt",
  PSA = "../02_harmonize_gwas/GCST90476277_PSA_hg38_qc.txt"
)

out_dir <- "./hyprcoloc_three_traits_four_loci"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

lead_ids <- c(
  KLK3 = "rs17632542",
  JAZF1 = "rs10486567",
  CLPTM1L = "rs414965",
  `NKX3-1` = "rs995432"
)

required_cols <- c(
  "SNP", "CHR", "BP", "effect_allele", "other_allele",
  "eaf", "beta", "se", "pval"
)

cat("Working directory:", getwd(), "\n")
cat("Start time:", format(Sys.time()), "\n")
cat("Reference panel:", REF, "\n")
cat("All results are generated from scratch, previous outputs are not used. \n\n")

stopifnot(file.exists(PLINK))
stopifnot(file.exists(paste0(REF, ".bed")))
stopifnot(file.exists(paste0(REF, ".bim")))
stopifnot(file.exists(paste0(REF, ".fam")))
stopifnot(all(file.exists(gwas_files)))

all_summaries <- list()

for (locus_name in names(lead_ids)) {

  cat("\n==================================================\n")
  cat("Processing locus:", locus_name, "\n")

  lead <- unname(lead_ids[[locus_name]])
  locus_prefix <- file.path(out_dir, locus_name)

  tryCatch({

    cat("Computing the LD variant set:", lead, "\n")

    ld_file <- paste0(locus_prefix, ".ld")
    plink_log <- paste0(locus_prefix, "_PLINK.log")

    status <- system2(
      PLINK,
      args = c(
        "--bfile", REF,
        "--ld-snp", lead,
        "--r2",
        "--ld-window-r2", "0.1",
        "--ld-window", "99999",
        "--ld-window-kb", "10000",
        "--out", locus_prefix
      ),
      stdout = plink_log,
      stderr = plink_log
    )

    if (status != 0 || !file.exists(ld_file)) {
      stop("PLINK run failed, check:", plink_log)
    }

    ld <- fread(ld_file, showProgress = FALSE)

    if (!"SNP_B" %in% names(ld)) {
      stop(
        "LD file Missing SNP_B columns; Actual columns Info:",
        paste(names(ld), collapse = ", ")
      )
    }

    ld_ids <- unique(c(lead, trimws(ld$SNP_B)))
    ld_ids <- ld_ids[!is.na(ld_ids) & ld_ids != ""]

    cat("LD Info rsID Info (Info lead):", length(ld_ids), "\n")

    fwrite(
      data.table(SNP = ld_ids),
      paste0(locus_prefix, "_ld_variant_rsids.csv")
    )

    gwas_regions <- list()

    for (trait in names(gwas_files)) {
      cat("Reading Info:", trait, "\n")

      d <- fread(
        gwas_files[[trait]],
        sep = "\t",
        header = TRUE,
        select = required_cols,
        showProgress = TRUE
      )

      d <- d[trimws(SNP) %in% ld_ids]

      d[, SNP := trimws(SNP)]
      d[, CHR := as.numeric(gsub("^chr", "", trimws(CHR)))]
      d[, BP := as.numeric(BP)]
      d[, effect_allele := toupper(trimws(effect_allele))]
      d[, other_allele := toupper(trimws(other_allele))]
      d[, eaf := as.numeric(eaf)]
      d[, beta := as.numeric(beta)]
      d[, se := as.numeric(se)]
      d[, pval := as.numeric(pval)]

      bad_numeric <- d[
        is.na(CHR) | is.na(BP) |
          is.na(beta) | is.na(se) | is.na(pval)
      ]

      if (nrow(bad_numeric) > 0) {
        cat(trait, "Failed numeric conversions: \n")
        print(bad_numeric)
        stop(trait, "Info Numeric conversion failed, records were not silently removed.")
      }

      d <- d[
        is.finite(beta) & is.finite(se) & is.finite(pval) &
          se > 0 & pval > 0 & pval <= 1
      ]

      d <- d[
        effect_allele %in% c("A", "C", "G", "T") &
          other_allele %in% c("A", "C", "G", "T")
      ]

      d[, palindromic :=
          (effect_allele == "A" & other_allele == "T") |
          (effect_allele == "T" & other_allele == "A") |
          (effect_allele == "C" & other_allele == "G") |
          (effect_allele == "G" & other_allele == "C")]

      cat(
        trait, ": records before QC =", nrow(d),
        "; palindromic records =", sum(d$palindromic), "\n"
      )

      d <- d[palindromic == FALSE]
      d[, palindromic := NULL]

      if (nrow(d) == 0) {
        stop(trait, "No variants remain after QC.")
      }

      gwas_regions[[trait]] <- copy(d)
      rm(d)
      gc()
    }

    for (trait in names(gwas_regions)) {
      if (!lead %in% gwas_regions[[trait]]$SNP) {
        stop(locus_name, "lead SNP not retained in", trait, "Info.")
      }
    }

    psa <- gwas_regions$PSA

    harm <- psa[, .(
      SNP,
      CHR,
      BP,
      EA = effect_allele,
      OA = other_allele,
      beta_PSA = beta,
      se_PSA = se,
      pval_PSA = pval
    )]

    harm[, variant_id := paste(
      SNP, CHR, BP, EA, OA,
      sep = ":"
    )]

    if (anyDuplicated(harm$variant_id)) {
      stop("PSA Info Duplicated variant_id.")
    }

    for (trait in c("BPH", "PCa")) {
      d <- gwas_regions[[trait]]

      m <- merge(
        harm[, .(SNP, variant_id, EA, OA)],
        d,
        by = "SNP",
        allow.cartesian = TRUE
      )

      same <- m$effect_allele == m$EA & m$other_allele == m$OA
      swapped <- m$effect_allele == m$OA & m$other_allele == m$EA
      keep_rows <- same | swapped

      m <- m[keep_rows]
      same <- same[keep_rows]

      if (nrow(m) == 0) {
        stop(trait, "and PSA No harmonizable allele records.")
      }

      m[, beta_aligned := ifelse(same, beta, -beta)]

      keep <- m[, .(
        variant_id,
        beta_aligned,
        se,
        pval
      )]

      setnames(
        keep,
        c("beta_aligned", "se", "pval"),
        c(
          paste0("beta_", trait),
          paste0("se_", trait),
          paste0("pval_", trait)
        )
      )

      if (anyDuplicated(keep$variant_id)) {
        stop(trait, "Duplicated after harmonization variant_id, records were not removed.")
      }

      harm <- merge(
        harm,
        keep,
        by = "variant_id",
        all = FALSE
      )
    }

    harm <- harm[
      is.finite(beta_BPH) & is.finite(se_BPH) & se_BPH > 0 &
        is.finite(beta_PCa) & is.finite(se_PCa) & se_PCa > 0 &
        is.finite(beta_PSA) & is.finite(se_PSA) & se_PSA > 0
    ]

    setorder(harm, CHR, BP, SNP)

    cat("Info Shared variants:", nrow(harm), "\n")

    if (nrow(harm) < 10) {
      stop("Fewer than 10, not running HyPrColoc.")
    }

    fwrite(
      harm,
      paste0(locus_prefix, "_harmonized_input.csv"),
      na = "NA"
    )

    effect_est <- as.matrix(
      harm[, .(
        BPH = beta_BPH,
        PCa = beta_PCa,
        PSA = beta_PSA
      )]
    )

    effect_se <- as.matrix(
      harm[, .(
        BPH = se_BPH,
        PCa = se_PCa,
        PSA = se_PSA
      )]
    )

    rownames(effect_est) <- harm$variant_id
    rownames(effect_se) <- harm$variant_id

    cat("Info rows HyPrColoc..\n")

    res <- hyprcoloc(
      effect.est = effect_est,
      effect.se = effect_se,
      binary.outcomes = c(TRUE, TRUE, TRUE),
      trait.names = c("BPH", "PCa", "PSA"),
      snp.id = harm$variant_id,
      ind.traits = TRUE
    )

    saveRDS(
      res,
      paste0(locus_prefix, "_hyprcoloc_full_results.rds")
    )

    result_table <- as.data.table(res$results)

    if (nrow(result_table) == 0) {
      stop("HyPrColoc returned no results.")
    }

    result_table[, `:=`(
      Locus = locus_name,
      Lead_SNP = lead,
      N_variants = nrow(harm)
    )]

    print(result_table)

    fwrite(
      result_table,
      paste0(locus_prefix, "_hyprcoloc_results.csv"),
      na = "NA"
    )

    if (!is.null(res$posteriors)) {
      fwrite(
        as.data.table(res$posteriors),
        paste0(locus_prefix, "_posterior_probabilities.csv"),
        na = "NA"
      )
    }

    if (!is.null(res$snpscores)) {
      snp_scores <- as.data.table(
        res$snpscores,
        keep.rownames = "variant_id"
      )

      fwrite(
        snp_scores,
        paste0(locus_prefix, "_per_variant_results.csv"),
        na = "NA"
      )
    }

    all_summaries[[locus_name]] <- copy(result_table)

    assoc <- harm[, .(
      marker = SNP,
      chr = as.integer(CHR),
      pos = as.integer(BP),
      z_1 = beta_BPH / se_BPH,
      z_2 = beta_PCa / se_PCa,
      z_3 = beta_PSA / se_PSA
    )]

    assoc <- assoc[
      is.finite(z_1) & is.finite(z_2) & is.finite(z_3)
    ]

    if (anyDuplicated(assoc$marker)) {
      dup_marker <- unique(
        assoc$marker[
          duplicated(assoc$marker) |
            duplicated(assoc$marker, fromLast = TRUE)
        ]
      )

      cat("Info marker: \n")
      print(dup_marker)
      stop(locus_name, "plotting input Duplicated marker, records were not removed.")
    }

    if (!lead %in% assoc$marker) {
      stop("lead SNP not in plotting input Info:", lead)
    }

    fwrite(
      assoc,
      paste0(locus_prefix, "_stack_assoc_plot_input.csv"),
      na = "NA"
    )

    cat("Plotting HyPrColoc stacked regional plot..\n")

    dat <- coloc_stack_assoc_plot(
      assoc = assoc,
      bfile_1000G = REF,
      SNP = lead,
      build = 38,
      traits = c("BPH", "PCa", "PSA"),
      save_plot = TRUE,
      width = 9,
      height = 10,
      plot_pdf = paste0(
        locus_name,
        "_HyPrColoc_stack_assoc"
      ),
      save_path = paste0(out_dir, "/")
    )

    saveRDS(
      dat,
      paste0(locus_prefix, "_stack_assoc_plot_object.rds")
    )

    cat("Completed locus:", locus_name, "\n")

  }, error = function(e) {
    cat("\n Info run failed:", locus_name, "\n")
    cat("Error:", conditionMessage(e), "\n")
  })
}

if (length(all_summaries) > 0) {
  summary_table <- rbindlist(
    all_summaries,
    use.names = TRUE,
    fill = TRUE
  )

  fwrite(
    summary_table,
    file.path(
      out_dir,
      "hyprcoloc_four_loci_main_summary.csv"
    ),
    na = "NA"
  )

  cat("\n========== Summary of completed loci ==========\n")
  print(summary_table)
}

cat("\n End time:", format(Sys.time()), "\n")
cat("Output directory:", normalizePath(out_dir, mustWork = FALSE), "\n")

sessionInfo()
