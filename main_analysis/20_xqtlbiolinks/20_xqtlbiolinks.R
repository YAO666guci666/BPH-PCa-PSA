library(data.table)
library(ggplot2)
library(xQTLbiolinks)

PLINK <- "/home/train/miniconda3/bin/plink"
PANEL <- "../reference/magma/g1000/g1000_eur"

if (!file.exists(PLINK)) stop("Not found plink:", PLINK)
if (!file.exists(paste0(PANEL, ".bim"))) {
  cat("PANEL Invalid prefix,", dirname(PANEL), "Actual contents: \n")
  print(list.files(dirname(PANEL)))
  stop("Set PANEL to one of the panel prefixes listed above")
}

get_LD <- function(target_snp, snp_vec, tag) {
  tmp_prefix <- file.path(tempdir(), paste0("ld_", gsub("[^A-Za-z0-9_]", "_", tag)))
  snp_file <- paste0(tmp_prefix, "_snps.txt")
  writeLines(unique(c(target_snp, snp_vec)), snp_file)
  cmd <- sprintf('"%s" --bfile "%s" --extract "%s" --r2 --ld-snp %s --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0 --out "%s"',
                 PLINK, PANEL, snp_file, target_snp, tmp_prefix)
  lg <- system(cmd, intern = TRUE)
  f <- paste0(tmp_prefix, ".ld")
  if (!file.exists(f)) {
    cat("  Warning:  ", target_snp, "LD could not be computed from the panel, plink log tail: \n")
    cat(paste(tail(lg, 4), collapse = "\n"), "\n")
    return(data.table(SNP_A = character(0), SNP_B = character(0), R2 = numeric(0)))
  }
  ld <- fread(f)
  snpLD <- rbind(ld[, .(SNP_A, SNP_B, R2)],
                 data.table(SNP_A = target_snp, SNP_B = target_snp, R2 = 1))
  cat("LD matched SNP:", nrow(snpLD) - 1, "(anchor:", target_snp, ")\n")
  snpLD
}

width_LocusZoom  <- 8
height_LocusZoom <- 6
wind_size        <- 50000
ref_pop          <- "g1000_eur"

GENES  <- c("RNASEH2B")
TISSUE <- "Prostate"

trait_map <- list(
  BPH = list(folder = "finngen_BPH",      gwas_fname = "finngen_BPH_hg38_qc.txt"),
  PSA = list(folder = "GCST90476277_PSA", gwas_fname = "GCST90476277_PSA_hg38_qc.txt"),
  PCa = list(folder = "ieu-b-4809_PC",    gwas_fname = "ieu-b-4809_PC_hg38_qc.txt")
)

OUT_DIR <- "./locus_5genes_prostate_Fast2TWAS"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TRAIT_FILE <- list(
  BPH = c(cav = "stable_eCAVIAR_finngen_BPH.sig.txt",       fen = "stable_finngen_BPH.fastenloc.sig.txt"),
  PSA = c(cav = "stable_eCAVIAR_GCST90476277_PSA.sig.txt",  fen = "stable_GCST90476277_PSA.fastenloc.sig.txt"),
  PCa = c(cav = "stable_eCAVIAR_ieu-b-4809_PC.sig.txt",      fen = "stable_ieu-b-4809_PC.fastenloc.sig.txt")
)

ef_file_path <- function(trait, method) {
  info  <- trait_map[[trait]]
  fname <- TRAIT_FILE[[trait]][if (method == "eCAVIAR") "cav" else "fen"]
  root  <- if (method == "eCAVIAR") "../10_process_ecaviar/" else "../13_process_fastenloc/"
  tried <- c(paste0(root, info$folder, "/", fname), paste0(root, fname))
  hit <- tried[file.exists(tried)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

pick_col <- function(nms, cands) { h <- cands[cands %in% nms]; if (length(h)) h[1] else NA_character_ }

read_ef_sig <- function(f, trait, method) {
  if (is.null(f) || !file.exists(f)) return(NULL)
  x <- fread(f)
  vcol <- if (method == "eCAVIAR") "CLPP" else "RCP"
  gcol <- pick_col(names(x), c("gene_symble", "gene_symbol", "gene"))
  lead <- pick_col(names(x), c("rsID_GWAS", "lead_snp", "lead_SNP", "rsID_lead"))
  rsrq <- pick_col(names(x), c("rsID_qtl", "rsID_QTL", "rsID", "QTL_rsid", "rsid_qtl"))
  qvid <- pick_col(names(x), c("qtl_Variant_id", "QTL_variant", "qtl_variant_id", "Variant_id", "QTL_Variant_id"))
  gvid <- pick_col(names(x), c("GWAS_variant", "gwas_variant", "GWAS_Variant_id"))
  need <- c(gcol, "Tissue", "QTL_Type", vcol)
  if (any(is.na(need)) || any(!need %in% names(x))) return(NULL)
  if (is.na(rsrq)) rsrq <- lead
  if (is.na(qvid)) qvid <- gvid
  if (is.na(lead)) return(NULL)
  x[, .(gkey = toupper(get(gcol)), trait, method,
        rsID_GWAS = get(lead), rsID_qtl = if (!is.na(rsrq)) get(rsrq) else get(lead),
        qtl_Variant_id = if (!is.na(qvid)) as.character(get(qvid)) else as.character(get(gvid)),
        Tissue, QTL_Type, val = as.numeric(get(vcol)))]
}

cat("Reading EF sig...\n")
ef_list <- list()
for (T in names(TRAIT_FILE)) {
  for (m in c("eCAVIAR", "fastenloc")) {
    f <- ef_file_path(T, m)
    if (!is.null(f)) {
      dt <- read_ef_sig(f, T, m)
      if (!is.null(dt)) ef_list[[length(ef_list)+1]] <- dt
    }
  }
}
ef <- rbindlist(ef_list, use.names = TRUE, fill = TRUE)
ef <- ef[gkey %in% toupper(GENES) & Tissue == TISSUE & !is.na(val)]
ef[, qt := ifelse(grepl("sQTL", QTL_Type, ignore.case = TRUE), "sQTL", "eQTL")]
ef_max <- ef[, .SD[which.max(val)], by = .(gkey, trait, method, Tissue, qt, rsID_qtl)]

cat("Total", nrow(ef_max), "records \n")
print(ef_max[, .(gkey, trait, method, Tissue, qt, rsID_GWAS, rsID_qtl, val = round(val, 4))])

pop <- toupper(gsub("g1000_", "", ref_pop))

for (i in seq_len(nrow(ef_max))) {

  row <- ef_max[i]
  gene     <- row$gkey
  trait    <- row$trait
  method   <- row$method
  snp      <- row$rsID_GWAS
  rsID     <- row$rsID_qtl
  qvid     <- row$qtl_Variant_id
  tissue   <- row$Tissue
  qtl_type <- row$QTL_Type
  val      <- row$val

  cat("\n[", i, "/", nrow(ef_max), "]", gene, "|", trait, "|", method, "| qtl:", rsID, "\n")

  snp_dir <- file.path(OUT_DIR, paste0(gene, "_", trait, "_", method))
  dir.create(snp_dir, showWarnings = FALSE, recursive = TRUE)

  info <- trait_map[[trait]]

  qtl_gwas_file <- paste0("../08_prepare_ecaviar/", info$folder, "/GWAS/",
                          snp, "_", info$folder, "_",
                          wind_size/1000, "kb_", ref_pop, "_markers_QTL_gwas_",
                          tissue, "_", qtl_type, ".all.txt")

  if (!file.exists(qtl_gwas_file)) {
    cat("Warning: QTL_GWAS does not exist:", basename(qtl_gwas_file), "\n"); next
  }

  qtl_gwas_data <- fread(qtl_gwas_file)
  if (!all(c("SNP", "pval_nominal", "pval") %in% names(qtl_gwas_data))) {
    cat("Warning: Missing columns \n"); next
  }

  ld_lead <- get_LD(snp, qtl_gwas_data$SNP, paste0(gene, "_", i, "_lead"))
  hl_lead <- snp
  if (nrow(ld_lead) == 0) {

    proxy <- qtl_gwas_data$SNP[which.min(qtl_gwas_data$pval)]
    cat("Warning: 03 LD anchor using within the window GWAS top SNP proxy:", proxy, "(instead of", snp, ")\n")
    ld_lead <- get_LD(proxy, qtl_gwas_data$SNP, paste0(gene, "_", i, "_leadproxy"))
    hl_lead <- proxy
  }
  if (identical(snp, rsID)) {
    ld_qtl <- ld_lead
  } else {
    ld_qtl <- get_LD(rsID, qtl_gwas_data$SNP, paste0(gene, "_", i, "_qtl"))
  }

  locusZoom_gwas_regiondata <- qtl_gwas_data[, .(SNP, pval)]
  locusZoom_qtl_regiondata  <- qtl_gwas_data[, .(SNP, pval_nominal)]

  xlab_text <- paste0(qtl_type, " –log10(P-value)")

  if (method == "eCAVIAR") {
    title1_base <- paste0("Lead GWAS variant: \n", snp)
    title2_base <- paste0("Variant with highest CLPP: \n", qvid, " (", rsID, ")")
    prefix2     <- paste0("max_CLPP_", rsID, "_", tissue)
    score_name  <- "CLPP"
  } else {
    title1_base <- paste0("Lead GWAS variant: \n", snp)
    title2_base <- paste0("Variant with highest RCP: \n", qvid, " (", rsID, ")")
    prefix2     <- paste0("max_RCP_", rsID, "_", tissue)
    score_name  <- "RCP"
  }

  title_p1 <- paste(gene, tissue, qtl_type, " \n", title1_base)
  tryCatch({
    p1 <- xQTLvisual_locusZoom(
      qtl_gwas_data[, c("SNP", "chr", "pos", "pval_nominal")],
      highlightSnp = hl_lead, population = pop, legend = TRUE, ylim = NULL,
      legend_position = "topleft", point_color = NULL,
      axis_text_size = 1.3, axis_title_size = 1.3, title_size = 1.4,
      ylab_text = xlab_text, title_text = title_p1,
      snpLD = ld_lead
    )
    print(p1)
    ggsave(paste0(snp_dir, "/", prefix2, "_", qtl_type, "_locusZoom_leadsnp.pdf"),
           p1, width = width_LocusZoom, height = height_LocusZoom)
    cat("   03\n")
  }, error = function(e) cat("   03:", e$message, "\n"))

  title_p2 <- paste(gene, tissue, qtl_type, " \n", title2_base)
  tryCatch({
    p2 <- xQTLvisual_locusZoom(
      qtl_gwas_data[, c("SNP", "chr", "pos", "pval_nominal")],
      highlightSnp = rsID, population = pop, legend = TRUE, ylim = NULL,
      legend_position = "topleft", point_color = NULL,
      axis_text_size = 1.3, axis_title_size = 1.3, title_size = 1.4,
      ylab_text = xlab_text, title_text = title_p2,
      snpLD = ld_qtl
    )
    print(p2)
    ggsave(paste0(snp_dir, "/", prefix2, "_", qtl_type, "_locusZoom_rsidqtl.pdf"),
           p2, width = width_LocusZoom, height = height_LocusZoom)
    cat("   04\n")
  }, error = function(e) cat("   04:", e$message, "\n"))

  title_p3 <- paste(gene, tissue, qtl_type, " vs\n", trait,
                    " GWAS (", method, " ", score_name, " = ", round(val, 4), ")\n",
                    qvid, ":", rsID)
  tryCatch({
    p3 <- xQTLvisual_locusCompare(
      locusZoom_qtl_regiondata, locusZoom_gwas_regiondata,
      highlightSnp = rsID, population = pop, legend = TRUE, legend_position = "topleft",
      axis_text_size = 1.3, axis_title_size = 1.3, title_size = 1.4,
      xlab_text = xlab_text, ylab_text = "GWAS –log10(P-value)", title_text = title_p3,
      snpLD = ld_qtl
    )
    print(p3)
    ggsave(paste0(snp_dir, "/", prefix2, "_", qtl_type, "_locusCompare.pdf"),
           p3, width = width_LocusZoom, height = height_LocusZoom)
    cat("   05\n")
  }, error = function(e) cat("   05:", e$message, "\n"))
}

cat("\n Script A Completed \n")

sessionInfo()
