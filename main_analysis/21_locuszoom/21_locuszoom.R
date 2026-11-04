library(data.table)
library(ggplot2)
library(Fast2TWAS)
library(qqman)
library(topr)

width_LocusZoom  <- 8
height_LocusZoom <- 6
wind_size        <- 300000
sig.thres        <- 5e-8
build            <- 38
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
  fname <- TRAIT_FILE[[trait]][if (method == "eCAVIAR") "cav" else "fen"]
  root  <- if (method == "eCAVIAR") "../10_process_ecaviar/" else "../13_process_fastenloc/"
  tried <- c(paste0(root, trait_map[[trait]]$folder, "/", fname),
             paste0(root, fname))
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

cat("Reading EF sig file...\n")
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

cat("Total", nrow(ef_max), "records to plot 01/02\n")
print(ef_max[, .(gkey, trait, method, Tissue, qt, rsID_GWAS, rsID_qtl, val = round(val, 4))])

get_snp_chrpos <- function(gwasfile, snp) {
  if (!file.exists(gwasfile)) return(c(CHR = NA_character_, BP = NA_character_))
  gwas <- fread(gwasfile)
  cn <- tolower(names(gwas))
  snp_col <- names(gwas)[which(cn %in% c("snp", "rsid", "rs", "markername"))[1]]
  chr_col <- names(gwas)[which(cn %in% c("chr", "chromosome", "chrom"))[1]]
  pos_col <- names(gwas)[which(cn %in% c("bp", "pos", "position"))[1]]
  if (any(is.na(c(snp_col, chr_col, pos_col)))) return(c(CHR = NA_character_, BP = NA_character_))
  m <- which(toupper(gwas[[snp_col]]) == toupper(snp))
  if (length(m) == 0) return(c(CHR = NA_character_, BP = NA_character_))
  c(CHR = as.character(gwas[[chr_col]][m[1]]), BP = as.character(gwas[[pos_col]][m[1]]))
}

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

  cat("\n[", i, "/", nrow(ef_max), "]", gene, "|", trait, "|", method,
      "| lead:", snp, "| qtl:", rsID, "\n")

  snp_dir <- file.path(OUT_DIR, paste0(gene, "_", trait, "_", method))
  dir.create(snp_dir, showWarnings = FALSE, recursive = TRUE)

  info <- trait_map[[trait]]

  gwas_file <- paste0("../02_harmonize_gwas/1/", info$gwas_fname)
  if (!file.exists(gwas_file)) { cat("Warning: GWAS does not exist:", gwas_file, "\n"); next }

  snp_pos <- get_snp_chrpos(gwas_file, snp)

  title1 <- paste0("Lead GWAS variant: \n", snp_pos["CHR"], ":", snp_pos["BP"], " (", snp, ")")
  if (is.na(snp_pos["CHR"])) title1 <- paste0("Lead GWAS variant: \n", snp)

  tryCatch({
    locusZoom_gassocplot2(
      gwasfile = gwas_file,
      wind_size = wind_size,
      SNP = snp,
      bfile = paste0("../reference/magma/g1000/", ref_pop),
      title = title1,
      sig.thres = sig.thres,
      width = width_LocusZoom,
      height = height_LocusZoom,
      build = build,
      outfile = paste0(snp_dir, "/leadSNPs_", trait, "_", snp)
    )
    cat("   01 leadSNPs_", trait, "_", snp, ".pdf\n")
  }, error = function(e) cat("   01:", e$message, "\n"))

  if (method == "eCAVIAR") {
    title2  <- paste0("Variant with highest CLPP: \n", qvid, " (", rsID, ")")
    prefix2 <- paste0("max_CLPP_", trait, "_", rsID, "_", tissue)
  } else {
    title2  <- paste0("Variant with highest RCP: \n", qvid, " (", rsID, ")\n(RCP = ", round(val, 3), ")")
    prefix2 <- paste0("max_RCP_", trait, "_", rsID, "_", tissue)
  }

  tryCatch({
    locusZoom_gassocplot2(
      gwasfile = gwas_file,
      wind_size = wind_size,
      SNP = rsID,
      bfile = paste0("../reference/magma/g1000/", ref_pop),
      title = title2,
      sig.thres = sig.thres,
      width = width_LocusZoom,
      height = height_LocusZoom,
      build = build,
      outfile = paste0(snp_dir, "/", prefix2)
    )
    cat("   02", prefix2, ".pdf\n")
  }, error = function(e) cat("   02:", e$message, "\n"))
}

cat("\n========== Script A Completed ==========\n")

sessionInfo()
