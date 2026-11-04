options(future.globals.maxSize = Inf)
if (!exists(".shared_pipeline_call")) rm(list = ls())
gc()
if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
library(data.table)
library(ieugwasr)
library(plinkbinr)
library(future.apply)
library(progressr)
plan(multisession, workers = max(1, floor(parallel::detectCores() * 0.5)))
if (!exists("inputFile")) inputFile <- c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")
wind_size <- 50000
ref_pop <- "g1000_eur"
trans <- FALSE
gtex_tissue_names <- fread("../reference/QTLEnrich/gtex_tissue_names.txt", header = FALSE)$V1
dir.create("logs", showWarnings = FALSE)
for (trait in inputFile) {
  GenomicRiskLoci_dir <- list.files(path = paste0("../07_fuma_gwas/", trait, "/"), pattern = "^FUMA", full.names = TRUE, include.dirs = TRUE)
  GenomicRiskLoci_dir <- GenomicRiskLoci_dir[file.info(GenomicRiskLoci_dir)$isdir][1]
  GenomicRiskLoci <- fread(paste0(GenomicRiskLoci_dir, "/GenomicRiskLoci.txt"))
  GenomicRiskLoci$chr <- gsub("^23$|^x$", "X", as.character(GenomicRiskLoci$chr))
  gwas_data <- fread(paste0("../02_harmonize_gwas/", trait, "_hg38_qc.txt"))
  gwas_data$CHR <- gsub("^23$|^x$", "X", as.character(gwas_data$CHR))
  trait_dir <- paste0(trait, "/")
  dir.create(trait_dir, showWarnings = FALSE, recursive = TRUE)
  trait_gwas_dir <- paste0(trait_dir, "GWAS")
  trait_corr_dir <- paste0(trait_dir, "corr")
  dir.create(trait_gwas_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(trait_corr_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(trait_dir, "eQTL"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(trait_dir, "sQTL"), showWarnings = FALSE, recursive = TRUE)
  progressr::with_progress({
    p <- progressor(along = GenomicRiskLoci$rsID)
    future_lapply(seq_along(GenomicRiskLoci$rsID), function(idx) {
      snp <- GenomicRiskLoci$rsID[idx]
      log_file <- sprintf("logs/%s_%s_%s.log", trait, snp, Sys.getpid())
      cat("", file = log_file)
      filter <- gwas_data[SNP == snp][1]
      if (nrow(filter) == 0) {
        p()
        return(NULL)
      }
      chrom <- filter$CHR
      pos <- filter$BP
      startbp <- pos - wind_size
      endbp <- pos + wind_size
      markers <- gwas_data[CHR == chrom & BP <= endbp & BP >= startbp]
      if (trans) {
        markers[, z_gwas := -beta / se]
        if ("eaf" %in% colnames(markers)) markers[, eaf := 1 - eaf]
        markers[, GWAS_variant := paste0("chr", CHR, "_", BP, "_", effect_allele, "_", other_allele, "_b38")]
      } else {
        markers[, z_gwas := beta / se]
        markers[, GWAS_variant := paste0("chr", CHR, "_", BP, "_", other_allele, "_", effect_allele, "_b38")]
      }
      for (qtl_type in c("eQTL", "sQTL")) {
        for (tissue in gtex_tissue_names) {
          qtl_file <- if (qtl_type == "eQTL") paste0("../reference/QTLEnrich/GTEx_Analysis_v8_eQTL/", tissue, ".v8.signif_variant_gene_pairs.txt.gz") else paste0("../reference/QTLEnrich/GTEx_Analysis_v8_sQTL/", tissue, ".v8.sqtl_signifpairs.txt.gz")
          qtl_dt <- fread(qtl_file)
          qtl_dt[, c("chr", "pos", "ref", "alt", "build") := tstrsplit(variant_id, "_")]
          markers_qtl <- qtl_dt[chr == paste0("chr", chrom) & as.integer(pos) <= endbp & as.integer(pos) >= startbp]
          setkey(markers_qtl, variant_id)
          markers_qtl <- unique(markers_qtl, by = "variant_id")
          markers_qtl_gwas <- merge(markers_qtl, markers, by.x = "variant_id", by.y = "GWAS_variant")
          if (nrow(markers_qtl_gwas) == 0) next
          corr <- try(ieugwasr::ld_matrix_local(markers_qtl_gwas$SNP, bfile = paste0("../reference/magma/g1000/", ref_pop), plink_bin = plinkbinr::get_plink_exe(), with_alleles = FALSE), silent = TRUE)
          if (inherits(corr, "try-error")) next
          final_SNP <- intersect(markers_qtl_gwas$SNP, colnames(corr))
          if (length(final_SNP) < 5) next
          markers_final <- markers_qtl_gwas[SNP %in% final_SNP]
          markers_final[, z_qtl := slope / slope_se]
          fwrite(markers_final, file = paste0(trait_gwas_dir, "/", snp, "_", trait, "_", wind_size / 1000, "kb_", ref_pop, "_markers_QTL_gwas_", tissue, "_", qtl_type, ".all.txt"), sep = "\t")
          fwrite(corr, file = paste0(trait_corr_dir, "/", snp, "_", trait, "_", wind_size / 1000, "kb_", ref_pop, "_markers_QTL_gwas_", tissue, "_", qtl_type, ".corr.txt"), sep = "\t", col.names = FALSE)
          fwrite(markers_final[, .(SNP, z_qtl)], file = paste0(trait_dir, qtl_type, "/", snp, "_", trait, "_", wind_size / 1000, "kb_", ref_pop, "_markers_QTL_gwas_", tissue, "_", qtl_type, ".qtl.txt"), sep = "\t", col.names = FALSE)
          fwrite(markers_final[, .(SNP, z_gwas)], file = paste0(trait_gwas_dir, "/", snp, "_", trait, "_", wind_size / 1000, "kb_", ref_pop, "_markers_QTL_gwas_", tissue, "_", qtl_type, ".gwas.txt"), sep = "\t", col.names = FALSE)
        }
      }
      p()
      NULL
    }, future.seed = TRUE)
  })
  gc()
}

sessionInfo()
