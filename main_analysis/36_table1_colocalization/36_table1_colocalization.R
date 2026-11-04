library(data.table)

pick_col <- function(nms, cands) { h <- cands[cands %in% nms]; if (length(h)) h[1] else NA_character_ }

parse_locus <- function(variant) {
  m <- regmatches(variant, regexpr("^chr[0-9XY]+_[0-9]+", variant))
  ifelse(nchar(m) > 0, sub("_", ":", m), NA_character_)
}

read_ef_sig <- function(f, trait, method) {
  x <- fread(f)
  vcol <- if (method == "eCAVIAR") "CLPP" else "RCP"
  gcol <- pick_col(names(x), c("gene_symble", "gene_symbol", "gene"))
  gid  <- pick_col(names(x), c("gene_id", "geneid", "GeneID", "ensembl_id", "phenotype_id"))
  lead <- pick_col(names(x), c("rsID_GWAS", "lead_snp", "lead_SNP", "rsID_lead"))
  rsrq <- pick_col(names(x), c("rsID_qtl", "rsID_QTL", "rsID", "QTL_rsid", "rsid_qtl"))
  qvid <- pick_col(names(x), c("qtl_Variant_id", "QTL_variant", "qtl_variant_id", "Variant_id", "QTL_Variant_id"))
  gvid <- pick_col(names(x), c("GWAS_variant", "gwas_variant", "GWAS_Variant_id"))
  hcol <- pick_col(names(x), c("hg38", "pos_hg38", "hg38_pos"))
  if (is.na(gcol) || is.na(lead)) stop("EF file Missing genes or lead SNP columns:\n", f)
  if (is.na(rsrq)) rsrq <- lead
  if (is.na(qvid)) qvid <- gvid
  if (is.na(hcol)) stop("EF file Missing hg38 columns, None Info lead GRCh38 Info:\n", f)
  out <- x[, .(Gene = get(gcol),
               gene_id = if (!is.na(gid)) as.character(get(gid)) else NA_character_,
               trait, method,
               lead = get(lead),
               rsID = if (!is.na(rsrq)) get(rsrq) else get(lead),
               variant = if (!is.na(qvid)) as.character(get(qvid)) else NA_character_,
               gwas_variant_raw = if (!is.na(gvid)) as.character(get(gvid)) else NA_character_,
               hg38 = as.numeric(get(hcol)),
               Tissue, QTL_Type, score = as.numeric(get(vcol)))]
  out[, locus := parse_locus(variant)]
  out[, chr := sub(":.*", "", locus)]
  out[, pos := suppressWarnings(as.numeric(sub(".*:", "", locus)))]
  out
}

ef_file_path <- function(trait, method) {
  trait_map <- list(
    BPH = list(folder = "finngen_BPH"),
    PSA = list(folder = "GCST90476277_PSA"),
    PCa = list(folder = "ieu-b-4809_PC")
  )
  TRAIT_FILE <- list(
    BPH = c(cav = "stable_eCAVIAR_finngen_BPH.sig.txt",       fen = "stable_finngen_BPH.fastenloc.sig.txt"),
    PSA = c(cav = "stable_eCAVIAR_GCST90476277_PSA.sig.txt",  fen = "stable_GCST90476277_PSA.fastenloc.sig.txt"),
    PCa = c(cav = "stable_eCAVIAR_ieu-b-4809_PC.sig.txt",     fen = "stable_ieu-b-4809_PC.fastenloc.sig.txt")
  )
  info  <- trait_map[[trait]]
  fname <- TRAIT_FILE[[trait]][if (method == "eCAVIAR") "cav" else "fen"]
  root  <- if (method == "eCAVIAR") "../10_process_ecaviar/" else "../13_process_fastenloc/"
  tried <- c(paste0(root, info$folder, "/", fname), paste0(root, fname))
  hit <- tried[file.exists(tried)]
  if (length(hit) == 0) stop("Not found EF file:\n", paste(tried, collapse = "\n  "))
  hit[1]
}

cat("Reading EF sig file...\n")
ef <- rbindlist(lapply(c("BPH", "PSA", "PCa"), function(T) {
  rbindlist(lapply(c("eCAVIAR", "fastenloc"), function(m) read_ef_sig(ef_file_path(T, m), T, m)))
}))
ef[, qt := ifelse(grepl("sQTL", QTL_Type, ignore.case = TRUE), "sQTL", "eQTL")]
cat("ef Info records:", nrow(ef), "\n\n")

best <- ef[, .SD[which.max(score)], by = .(trait, lead, Gene, Tissue, qt, method)]
best_cav <- best[method == "eCAVIAR", .(trait, lead, Gene, Tissue, qt,
                                        CLPP = score, gid_cav = gene_id, rsID_cav = rsID, variant_cav = variant,
                                        locus_cav = locus, chr, pos,
                                        gwas_variant_raw = gwas_variant_raw[1])]
best_fen <- best[method == "fastenloc", .(trait, lead, Gene, Tissue, qt,
                                          RCP = score, gid_fen = gene_id, rsID_fen = rsID, variant_fen = variant,
                                          locus_fen = locus, chr_fen = chr, pos_fen = pos,
                                          gwas_variant_raw = gwas_variant_raw[1])]
wide_new <- merge(best_cav, best_fen,
                  by = c("trait", "lead", "Gene", "Tissue", "qt", "gwas_variant_raw"), all = TRUE)

wide_new[, chr := fcoalesce(chr, chr_fen)]
wide_new[, pos := fcoalesce(pos, pos_fen)]
wide_new[, c("chr_fen", "pos_fen") := NULL]

wide_new[, is_prostate := grepl("prostate", Tissue, ignore.case = TRUE)]
wide_new[, single_method := (!is.na(CLPP) & CLPP >= 0.01) | (!is.na(RCP) & RCP >= 0.1)]

MAX_GTEX_GENES <- 3
WIND_SIZE <- 50000

GTF_FILE <- "gencode.v44.annotation.gtf.gz"
if (!file.exists(GTF_FILE)) {
  cat("Working directory None", GTF_FILE, ", Info Info (GRCh38, GENCODE v44)...\n")
  download.file("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz",
                destfile = GTF_FILE, mode = "wb")
}
if (!file.exists(GTF_FILE)) stop("GTF Error: Info Info gencode.v44.annotation.gtf.gz Info Working directory")

gtf <- fread(GTF_FILE, sep = "\t", header = FALSE, comment.char = "#",
             col.names = c("seqname", "source", "feature", "start", "end",
                           "score", "strand", "frame", "attribute"))
gtf <- gtf[feature == "gene"]
gene_coord <- gtf[, .(
  gene_id   = sub("\\..*$", "", sub(".*gene_id \"([^\"]+)\".*", "\\1", attribute)),
  gene_name = sub(".*gene_name \"([^\"]+)\".*", "\\1", attribute),
  chr = sub("^chr", "", seqname),
  start = start,
  end   = end
)]
gene_coord <- gene_coord[chr %in% c(as.character(1:22), "X", "Y", "M")]
gene_coord <- unique(gene_coord, by = "gene_id")
gene_coord <- gene_coord[!grepl("^ENSG", gene_name)]
cat("GENCODE v44 genes Info:", nrow(gene_coord), "genes \n")

make_gene_tag <- function(dt, max_genes = 3) {
  if (nrow(dt) == 0) return("-")
  dt[, display_score := pmax(CLPP / 0.01, RCP / 0.1, na.rm = TRUE)]
  dt[is.infinite(display_score), display_score := NA_real_]
  gene_info <- dt[, .(display_score = max(display_score, na.rm = TRUE)), by = Gene]
  gene_info <- gene_info[order(-display_score)]
  genes <- gene_info$Gene
  n <- length(genes)
  parts <- sapply(genes[seq_len(min(n, max_genes))], function(g) {
    sub <- dt[dt$Gene == g, ]
    qt_chars <- sort(unique(substr(sub$qt, 1, 1)))
    qt_str <- paste(qt_chars, collapse = ",")
    paste0(g, " (", qt_str, ")")
  })
  out <- paste(parts, collapse = "; ")
  if (n > max_genes) out <- paste0(out, "; +", n - max_genes, " more")
  out
}

OUT <- "./Table1_colocalization"
dir.create(OUT, showWarnings = FALSE)

lead_snps <- ef[, .(chr = chr[1],
                    pos = hg38[!is.na(hg38)][1],
                    gwas_variant_raw = gwas_variant_raw[1]),
                by = .(trait, lead)]
lead_snps <- lead_snps[!is.na(pos)]
setorder(lead_snps, trait, chr, pos)
cat("Lead SNP Info:", nrow(lead_snps), "\n")

make_hamel_by_lead <- function(anchor) {
  leads <- lead_snps[trait == anchor]
  if (nrow(leads) == 0) return(NULL)

  rbindlist(lapply(seq_len(nrow(leads)), function(i) {
    r <- leads[i]
    ld <- r$lead
    chr0 <- r$chr
    pos0 <- r$pos
    pmin0 <- pos0 - WIND_SIZE
    pmax0 <- pos0 + WIND_SIZE

    genes_cav <- ef[trait == anchor & lead == ld & method == "eCAVIAR" & score >= 0.01, unique(Gene)]
    genes_fen <- ef[trait == anchor & lead == ld & method == "fastenloc" & score >= 0.1, unique(Gene)]
    n_total <- length(union(genes_cav, genes_fen))

    dm_prostate <- wide_new[trait == anchor & lead == ld &
                              is_prostate == TRUE & single_method == TRUE]
    prostate_tag <- make_gene_tag(dm_prostate, max_genes = MAX_GTEX_GENES)

    dm_hc <- wide_new[trait == anchor & lead == ld & single_method == TRUE]
    tag <- make_gene_tag(dm_hc, max_genes = MAX_GTEX_GENES)

    g_sub <- gene_coord[chr == sub("^chr", "", chr0)]
    if (nrow(g_sub) > 0) {
      g_sub[, dist := pmax(start - pos0, pos0 - end, 0)]
      nearest <- g_sub[order(dist, gene_name)][1, gene_name]
    } else nearest <- ""

    gv <- sub("_b37$|_b38$|_GRCh37$|_GRCh38$|_hg19$|_hg38$", "", r$gwas_variant_raw)
    if (grepl("^chr[0-9XY]+[:_][0-9]+[:_][ACGT]+[:_][ACGT]+$", gv)) {
      alleles <- sub("^chr[0-9XY]+[:_][0-9]+[:_]", "", gv)
      lead_fmt <- paste0(chr0, ":", pos0, ":", gsub("_", ":", alleles))
    } else if (grepl("^rs\\d+$", gv)) {
      lead_fmt <- gv
    } else {
      lead_fmt <- paste0(chr0, ":", pos0)
    }

    out <- data.table(
      Lead_GWAS_variant = lead_fmt,
      RS_ID = ld,
      Nearest_Gene = nearest,
      Total_colocalizing_esGenes = n_total,
      Prostate_esQTLs = prostate_tag,
      Colocalizing_GTEx_esQTLs = tag
    )

    for (pt in setdiff(c("PSA", "BPH", "PCa"), anchor)) {
      has_cross <- nrow(L_cross_pool[trait == pt & chr == chr0 &
                                       pos >= pmin0 & pos <= pmax0]) > 0
      out[[paste0("Colocalizing_with_", pt)]] <- if (has_cross) "Yes" else "No"
    }

    out
  }))
}

hamel_lead <- setNames(lapply(c("PSA", "BPH", "PCa"), make_hamel_by_lead),
                       c("PSA", "BPH", "PCa"))

tb <- hamel_lead[["PSA"]]

cat("\n=== PSA Hamel Info (v12_final_bom_nearest4) ===\n")
cat("Info lead Info:", nrow(tb), "\n")
cat("Info Colocalizing_GTEx_esQTLs Info:", nrow(tb[Colocalizing_GTEx_esQTLs != "-"]), "\n")
cat("Colocalizing_with_PCa = Yes:", nrow(tb[Colocalizing_with_PCa == "Yes"]), "\n")
cat("Colocalizing_with_BPH = Yes:", nrow(tb[Colocalizing_with_BPH == "Yes"]), "\n")

main <- tb[Colocalizing_with_PCa == "Yes" | Colocalizing_with_BPH == "Yes"]
setorder(main, -Total_colocalizing_esGenes)

cat("\n Main-text Info (PCa=' Info ' | BPH=' Info '):", nrow(main), "Info \n")

print(main[, .(Lead_GWAS_variant, RS_ID, Nearest_Gene,
               Total_colocalizing_esGenes, Prostate_esQTLs,
               Colocalizing_GTEx_esQTLs, Colocalizing_with_PCa, Colocalizing_with_BPH)])

fwrite(tb, file.path(OUT, "Supplementary_PSA_Hamel_byLead_all_loci.csv"), bom = TRUE)
fwrite(main, file.path(OUT, "Table1_PSA_Hamel_byLead_main.csv"), bom = TRUE)

fwrite(hamel_lead[["BPH"]], file.path(OUT, "Supplementary_BPH_Hamel_byLead.csv"), bom = TRUE)
fwrite(hamel_lead[["PCa"]], file.path(OUT, "Supplementary_PCa_Hamel_byLead.csv"), bom = TRUE)

cat("\n Output directory:", normalizePath(OUT), "\n")

main_en <- main[, .(Lead_GWAS_variant, RS_ID, Nearest_Gene,
                    Total_colocalizing_esGenes, Prostate_esQTLs,
                    Colocalizing_with_PCa, Colocalizing_with_BPH)]

setnames(main_en,
         old = c("Lead_GWAS_variant", "RS_ID", "Nearest_Gene",
                 "Total_colocalizing_esGenes", "Prostate_esQTLs",
                 "Colocalizing_with_PCa", "Colocalizing_with_BPH"),
         new = c("Lead GWAS variant", "Lead SNP", "Nearest gene",
                 "Total colocalizing e/sGenes", "Prostate e/sQTLs",
                 "Colocalizing with PCa", "Colocalizing with BPH"))

fwrite(main_en, file.path(OUT, "Table1_PSA_Hamel_English.csv"), bom = TRUE)

footnote_text <- "Table 1. Colocalizing e/sGenes at PSA GWAS loci with shared genetic architecture for prostate cancer and benign prostatic hyperplasia.

Column descriptions:
  Lead GWAS variant          Chromosomal position of the lead SNP (chr:pos:ref:alt, GRCh38).
  Lead SNP                   rsID of the lead GWAS variant.
  Nearest gene               Gene with an approved symbol physically closest to the lead SNP
                             (GENCODE v44 gene annotation, GRCh38; distance measured to the
                             gene body, 0 if the variant falls within the gene; genes without
                             an approved symbol, i.e. ENSG identifiers in GENCODE v44, were
                             excluded).
  Total colocalizing e/sGenes
                             Number of genes with colocalization evidence (eCAVIAR CLPP >= 0.01
                             or fastENLOC RCP >= 0.1) at the locus.
  Prostate e/sQTLs           Genes showing colocalization in prostate tissue (e, eQTL; s, sQTL),
                             i.e. genes with eCAVIAR CLPP >= 0.01 or fastENLOC RCP >= 0.1 in
                             prostate tissue at the same lead SNP.
  Colocalizing with PCa      Cross-trait colocalization with prostate cancer GWAS. Yes = genes
                             with eCAVIAR CLPP >= 0.01 or fastENLOC RCP >= 0.1 present within
                             +/- 50 kb of the lead SNP; No = no colocalizing gene in the region.
  Colocalizing with BPH      Cross-trait colocalization with benign prostatic hyperplasia GWAS.
                             Yes = same criterion as for PCa; No = no colocalizing gene.

Methods summary:
  Colocalization:       eCAVIAR CLPP >= 0.01 or fastENLOC RCP >= 0.1.
  Cross-trait matching: Lead SNP +/- 50 kb window; colocalization status determined by
                        presence of colocalizing genes in the matched region.
"

writeLines(footnote_text, con = file.path(OUT, "Table1_footnote.txt"))
cat("Info Saved \n")

sessionInfo()
