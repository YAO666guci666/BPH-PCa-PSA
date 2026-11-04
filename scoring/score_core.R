required_config <- c("OUT_DIR", "TOP_N", "CLPP_CUT", "RCP_CUT", "EF_FILES",
                     "SCDRS_FILES", "PAIRS", "CAPTION_TEXT", "EXPECTED_MAIN_GENES", "O4G_FILE")
missing_config <- required_config[!vapply(required_config, exists, logical(1), inherits = TRUE)]
if (length(missing_config) > 0) stop("Missing Info:", paste(missing_config, collapse = ", "))
library(data.table)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

read_ef <- function(T) {
  fs <- EF_FILES[[T]]
  res <- list()
  for (k in names(fs)) {
    f <- fs[k]
    if (!file.exists(f)) stop("Not found EF file:", f)
    x <- fread(f)
    gcol <- intersect(c("gene_symble", "gene_symbol", "gene"), names(x))[1]
    vcol <- if (k == "cav") "CLPP" else "RCP"
    need <- c(gcol, "Tissue", "QTL_Type", vcol)
    miss <- setdiff(need, names(x))
    if (length(miss) > 0) stop(f, "Info columns:", paste(miss, collapse = ", "))
    x <- x[, .(gkey = toupper(get(gcol)), Tissue, QTL_Type,
               val = suppressWarnings(as.numeric(get(vcol))))]
    x <- x[!is.na(gkey) & gkey != ""]
    x[, is_pro := grepl("prostate", Tissue, ignore.case = TRUE)]
    res[[paste0(k, "_cross")]] <- x[is_pro == FALSE & !is.na(val),
                                    .(cross_n = uniqueN(Tissue), vmax = max(val)), by = gkey]
    res[[paste0(k, "_pro")]] <- x[is_pro == TRUE & !is.na(val),
                                  .(pval = max(val)), by = gkey]
  }
  all_g <- unique(unlist(lapply(res, function(z) z$gkey)))
  cat("  EF", T, "genes Info:", length(all_g), "\n")
  c(res, list(all = all_g))
}

read_o4g <- function() {
  f <- O4G_FILE
  if (!file.exists(f)) stop("Not found Open4Gene results file:", f)
  d <- if (grepl("[.]xlsx$", f, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Reading xlsx Info readxl Info")
    as.data.table(readxl::read_excel(f))
  } else fread(f)
  gcol <- intersect(c("Gene","gene","SYMBOL","symbol","GENE","gene_symbol"), names(d))[1]
  if (is.na(gcol)) stop(f, "Info gene columns, Actual columns:", paste(names(d), collapse = " | "))
  g <- unique(toupper(d[[gcol]])); g <- g[!is.na(g) & g != ""]
  cat("  Open4Gene <-", f, "genes Info:", length(g), "\n")
  g
}

read_placo <- function(pair) {
  base_pair <- sub("_6319$", "", pair)
  candidates <- unique(c(paste0(pair, "_genes.xlsx"), paste0(pair, "_genes.csv"),
                  paste0(toupper(pair), "_genes.xlsx"), paste0(toupper(pair), "_genes.csv"),
                  paste0(base_pair, "_genes.xlsx"), paste0(base_pair, "_genes.csv"),
                  paste0(toupper(base_pair), "_genes.xlsx"), paste0(toupper(base_pair), "_genes.csv")))
  for (f in candidates) {
    if (!file.exists(f)) next
    d <- if (grepl("[.]csv$", f, ignore.case = TRUE)) tryCatch(fread(f), error = function(e) NULL)
    else if (requireNamespace("readxl", quietly = TRUE)) tryCatch(as.data.table(readxl::read_excel(f)), error = function(e) NULL)
    else NULL
    if (!is.null(d) && nrow(d) > 0) {
      gcol <- intersect(c("Gene","gene","SYMBOL","symbol","GENE"), names(d))[1]
      if (is.na(gcol)) next
      g <- unique(toupper(d[[gcol]])); g <- g[!is.na(g) & g != ""]
      cat("  PLACO", pair, "<-", f, "genes Info:", length(g), "\n")
      return(g)
    }
  }
  stop("Not found PLACO Info file:", pair)
}

read_scdrs <- function(pair, top_n = TOP_N) {
  f <- unname(SCDRS_FILES[pair])
  if (!file.exists(f)) stop("Not found scDRS file:", f)
  cor_csv <- sub("_top1000_genes[.]txt$", "_genomewide_cor.csv", f)
  if (file.exists(cor_csv)) {
    d <- fread(cor_csv)
    gcol <- intersect(c("gene","symbol","gene_name","Gene","SYMBOL","GENE"), names(d))[1]
    if (is.na(gcol) || !"cor" %in% names(d)) stop(cor_csv, "Info gene columns and cor columns, Actual columns:", paste(names(d), collapse = " | "))
    d[, gkey := toupper(get(gcol))]
    d[, cor := suppressWarnings(as.numeric(cor))]
    d <- d[!is.na(gkey) & gkey != "" & !is.na(cor)]
    setorder(d, -cor); d <- d[!duplicated(gkey)]
    top <- if (nrow(d) > top_n) d[1:top_n] else d
    cat("  scDRS", pair, "<-", cor_csv, " top", top_n, "threshold cor =", round(min(top$cor), 4), "\n")
    return(list(hit = top$gkey, cor_all = d[, .(gkey, cor)]))
  }
  d <- fread(f, header = FALSE)
  if (ncol(d) != 1) stop(f, "Info columns genes Info (None header)txt or genomewide_cor.csv, Actual columns Info:", ncol(d))
  g <- unique(toupper(d[[1]])); g <- g[!is.na(g) & g != ""]
  cat("  scDRS", pair, "<-", f, " top", length(g), "(no correlation CSV, cor set NA, shown as -)\n")
  list(hit = g, cor_all = data.table(gkey = g, cor = NA_real_))
}

master <- list()
o4g_global <- read_o4g()

for (cfg in PAIRS) {
  cat("\n==========", cfg$pair, "==========\n")
  ef1 <- read_ef(cfg$T1); ef2 <- read_ef(cfg$T2)
  pl  <- read_placo(cfg$pair)
  sc  <- read_scdrs(cfg$pair)

  genes <- sort(unique(c(ef1$all, ef2$all, o4g_global, pl, sc$hit)))
  genes <- genes[!is.na(genes) & genes != ""]
  cat("All genes =", length(genes), "\n")
  tab <- data.table(Gene = genes)

  tab[, cross_n_cav1 := ef1$cav_cross$cross_n[match(tab$Gene, ef1$cav_cross$gkey)]]
  tab[, cross_n_fen1 := ef1$fen_cross$cross_n[match(tab$Gene, ef1$fen_cross$gkey)]]
  tab[, cross_n_tissues_T1 := pmax(cross_n_cav1, cross_n_fen1, na.rm = TRUE)]
  tab[is.infinite(cross_n_tissues_T1), cross_n_tissues_T1 := NA]
  tab[, max_CLPP_T1 := ef1$cav_cross$vmax[match(tab$Gene, ef1$cav_cross$gkey)]]
  tab[, max_RCP_T1  := ef1$fen_cross$vmax[match(tab$Gene, ef1$fen_cross$gkey)]]
  tab[, cross_n_cav2 := ef2$cav_cross$cross_n[match(tab$Gene, ef2$cav_cross$gkey)]]
  tab[, cross_n_fen2 := ef2$fen_cross$cross_n[match(tab$Gene, ef2$fen_cross$gkey)]]
  tab[, cross_n_tissues_T2 := pmax(cross_n_cav2, cross_n_fen2, na.rm = TRUE)]
  tab[is.infinite(cross_n_tissues_T2), cross_n_tissues_T2 := NA]
  tab[, max_CLPP_T2 := ef2$cav_cross$vmax[match(tab$Gene, ef2$cav_cross$gkey)]]
  tab[, max_RCP_T2  := ef2$fen_cross$vmax[match(tab$Gene, ef2$fen_cross$gkey)]]

  tab[, prostate_CLPP_T1 := ef1$cav_pro$pval[match(tab$Gene, ef1$cav_pro$gkey)]]
  tab[, prostate_RCP_T1  := ef1$fen_pro$pval[match(tab$Gene, ef1$fen_pro$gkey)]]
  tab[, prostate_CLPP_T2 := ef2$cav_pro$pval[match(tab$Gene, ef2$cav_pro$gkey)]]
  tab[, prostate_RCP_T2  := ef2$fen_pro$pval[match(tab$Gene, ef2$fen_pro$gkey)]]

  tab[, cav1 := as.integer(!is.na(max_CLPP_T1) & max_CLPP_T1 >= CLPP_CUT)]
  tab[, cav2 := as.integer(!is.na(max_CLPP_T2) & max_CLPP_T2 >= CLPP_CUT)]
  tab[, fen1 := as.integer(!is.na(max_RCP_T1)  & max_RCP_T1  >= RCP_CUT)]
  tab[, fen2 := as.integer(!is.na(max_RCP_T2)  & max_RCP_T2  >= RCP_CUT)]
  tab[, G_ef := cav1 + cav2 + fen1 + fen2]
  tab[, Open4Gene := as.integer(tab$Gene %in% o4g_global)]
  tab[, G_score := G_ef + Open4Gene]

  tab[, prostate_T1 := as.integer(
    (!is.na(prostate_CLPP_T1) & prostate_CLPP_T1 >= CLPP_CUT) |
      (!is.na(prostate_RCP_T1)  & prostate_RCP_T1  >= RCP_CUT))]
  tab[, prostate_T2 := as.integer(
    (!is.na(prostate_CLPP_T2) & prostate_CLPP_T2 >= CLPP_CUT) |
      (!is.na(prostate_RCP_T2)  & prostate_RCP_T2  >= RCP_CUT))]
  tab[, scDRS := as.integer(tab$Gene %in% sc$hit)]
  tab[, scDRS_cor := sc$cor_all$cor[match(tab$Gene, sc$cor_all$gkey)]]
  tab[, T_score := prostate_T1 + prostate_T2 + scDRS]

  tab[, PLACO_MAGMA := as.integer(tab$Gene %in% pl)]

  tab[, Total := T_score + G_score]
  tmax <- max(tab$Total[tab$T_score >= 2L])
  tab[, Candidate := T_score >= 2L & Total == tmax]
  tab[, Route := fcase(
    T_score >= 2L & Total == tmax, "Maximum Total score in pair",
    default = "-")]

  cat("Candidates:", sum(tab$Candidate), ", Maximum Total =", tmax, ": ",
      paste(tab$Gene[tab$Candidate], collapse = ", "), "\n")

  keep_cols <- c("Gene", "PLACO_MAGMA",
                 "cav1", "cav2", "fen1", "fen2",
                 "Open4Gene",
                 "prostate_T1", "prostate_T2", "scDRS",
                 "Total", "Candidate", "Route",
                 "cross_n_tissues_T1", "cross_n_tissues_T2",
                 "max_CLPP_T1", "max_RCP_T1", "max_CLPP_T2", "max_RCP_T2",
                 "prostate_CLPP_T1", "prostate_RCP_T1", "prostate_CLPP_T2", "prostate_RCP_T2",
                 "scDRS_cor")
  d_out <- tab[, keep_cols, with = FALSE]
  setorder(d_out, -Candidate, -Total, Gene)
  fwrite(d_out, file.path(OUT_DIR, paste0("Table_", cfg$pair, "_unified_scores.csv")))
  fwrite(d_out[Candidate == TRUE], file.path(OUT_DIR, paste0("Table_", cfg$pair, "_candidate_genes.csv")))
  cat("Written:", cfg$pair, "\n")

  fmt_pair <- function(v1, v2) {
    s1 <- ifelse(is.na(v1), "-", sprintf("%.3f", v1))
    s2 <- ifelse(is.na(v2), "-", sprintf("%.3f", v2))
    out <- paste0(s1, " / ", s2)
    out[s1 == "-" & s2 == "-"] <- "-"
    out
  }
  d_cand <- tab[Candidate == TRUE]
  T1n <- cfg$T1
  main_tab <- data.frame(
    Gene = d_cand$Gene,
    "Colocalization_T1_non_prostate_CLPP_RCP" = fmt_pair(d_cand$max_CLPP_T1, d_cand$max_RCP_T1),
    "Colocalization_PSA_non_prostate_CLPP_RCP" = fmt_pair(d_cand$max_CLPP_T2, d_cand$max_RCP_T2),
    "Colocalization_T1_prostate_CLPP_RCP" = fmt_pair(d_cand$prostate_CLPP_T1, d_cand$prostate_RCP_T1),
    "Colocalization_PSA_prostate_CLPP_RCP" = fmt_pair(d_cand$prostate_CLPP_T2, d_cand$prostate_RCP_T2),
    "Open4Gene" = ifelse(d_cand$Open4Gene == 1, "Y", "-"),
    "scDRS_correlation" = ifelse(d_cand$scDRS == 1 & !is.na(d_cand$scDRS_cor),
                                 sprintf("%.3f", d_cand$scDRS_cor), "-"),
    "Total_Score" = d_cand$Total,
    check.names = FALSE
  )
  names(main_tab) <- gsub("_T1", paste0("_", T1n), names(main_tab), fixed = TRUE)
  main_tab <- main_tab[order(-main_tab[["Total_Score"]]), ]
  fwrite(main_tab, file.path(OUT_DIR, paste0("Table2_", cfg$pair, "_main_table.csv")))
  cat("Main table:", nrow(main_tab), "rows ->", paste0("Table2_", cfg$pair, "_main_table.csv"), "\n")

  master[[cfg$pair]] <- d_out[, pair := cfg$pair]
}

fwrite(rbindlist(master, fill = TRUE), file.path(OUT_DIR, "Table_master_unified_scores.csv"))

cat("\n========== Building Main-text Table 2(combined)==========\n")

bph_psa <- master[[PAIRS[[1]]$pair]]
pca_psa <- master[[PAIRS[[2]]$pair]]
setkey(bph_psa, Gene); setkey(pca_psa, Gene)

all_genes <- sort(unique(c(
  bph_psa[Total >= 7, Gene],
  pca_psa[Total >= 7, Gene]
)))
cat("Main-text Table 2 genes Info:", length(all_genes), "\n")
if (!is.na(EXPECTED_MAIN_GENES) && length(all_genes) != EXPECTED_MAIN_GENES) warning("main analysis Info genes Info", EXPECTED_MAIN_GENES, ", Check input version, Info genes")

safe_get <- function(dt, col) {
  if (nrow(dt) == 0) return(NA_real_)
  val <- dt[[col]]
  if (length(val) == 0) return(NA_real_)
  return(val[1])
}

fmt_triple <- function(bph, psa, pca) {
  b <- ifelse(is.na(bph), "–", sprintf("%.3f", bph))
  p <- ifelse(is.na(psa), "–", sprintf("%.3f", psa))
  c <- ifelse(is.na(pca), "–", sprintf("%.3f", pca))
  if (b == "–" && p == "–" && c == "–") return("–")
  paste(b, p, c, sep = " / ")
}

tab2 <- data.table(Gene = all_genes)

for (i in seq_len(nrow(tab2))) {
  g <- tab2$Gene[i]
  b1 <- bph_psa[g]
  p1 <- pca_psa[g]

  bph_cav <- safe_get(b1, "max_CLPP_T1")
  psa_cav <- ifelse(!is.na(safe_get(b1, "max_CLPP_T2")), safe_get(b1, "max_CLPP_T2"), safe_get(p1, "max_CLPP_T2"))
  pca_cav <- safe_get(p1, "max_CLPP_T1")

  bph_fen <- safe_get(b1, "max_RCP_T1")
  psa_fen <- ifelse(!is.na(safe_get(b1, "max_RCP_T2")), safe_get(b1, "max_RCP_T2"), safe_get(p1, "max_RCP_T2"))
  pca_fen <- safe_get(p1, "max_RCP_T1")

  o4g_bph <- safe_get(b1, "Open4Gene")
  o4g_pca <- safe_get(p1, "Open4Gene")
  o4g <- ifelse((!is.na(o4g_bph) && o4g_bph == 1) || (!is.na(o4g_pca) && o4g_pca == 1), 1, 0)

  sc_bph_hit <- safe_get(b1, "scDRS")
  sc_pca_hit <- safe_get(p1, "scDRS")
  sc_bph_cor <- safe_get(b1, "scDRS_cor")
  sc_pca_cor <- safe_get(p1, "scDRS_cor")

  bph_total <- safe_get(b1, "Total")
  pca_total <- safe_get(p1, "Total")

  tab2[i, `:=`(
    eCAVIAR_Non_Prostate_Maximum_CLPP = fmt_triple(bph_cav, psa_cav, pca_cav),
    fastENLOC_Non_Prostate_Maximum_RCP = fmt_triple(bph_fen, psa_fen, pca_fen),
    Open4Gene = ifelse(o4g == 1, "✓", "–"),
    scDRS_BPH_PSA_Correlation = ifelse(is.na(sc_bph_hit) || sc_bph_hit == 0 || is.na(sc_bph_cor), "–", sprintf("%.3f", sc_bph_cor)),
    scDRS_PCa_PSA_Correlation = ifelse(is.na(sc_pca_hit) || sc_pca_hit == 0 || is.na(sc_pca_cor), "–", sprintf("%.3f", sc_pca_cor)),
    BPH_PSA_Total_Score = ifelse(is.na(bph_total), "–", as.character(bph_total)),
    PCa_PSA_Total_Score = ifelse(is.na(pca_total), "–", as.character(pca_total))
  )]
}

tab2[, BPH_Num := ifelse(BPH_PSA_Total_Score == "–", 0, as.numeric(BPH_PSA_Total_Score))]
tab2[, PCa_Num := ifelse(PCa_PSA_Total_Score == "–", 0, as.numeric(PCa_PSA_Total_Score))]
setorder(tab2, -PCa_Num, -BPH_Num, Gene)
tab2[, c("BPH_Num", "PCa_Num") := NULL]

sc1_col <- paste0("scDRS_", PAIRS[[1]]$pair, "_Correlation")
sc2_col <- paste0("scDRS_", PAIRS[[2]]$pair, "_Correlation")
total1_col <- paste0(PAIRS[[1]]$pair, "_Total_Score")
total2_col <- paste0(PAIRS[[2]]$pair, "_Total_Score")
setnames(tab2,
         c("scDRS_BPH_PSA_Correlation", "scDRS_PCa_PSA_Correlation",
           "BPH_PSA_Total_Score", "PCa_PSA_Total_Score"),
         c(sc1_col, sc2_col, total1_col, total2_col))
setcolorder(tab2, c("Gene",
                    "eCAVIAR_Non_Prostate_Maximum_CLPP",
                    "fastENLOC_Non_Prostate_Maximum_RCP",
                    "Open4Gene", sc1_col, sc2_col, total1_col, total2_col))

fwrite(tab2, file.path(OUT_DIR, "Table2_Main_Text.csv"))
cat("Main-text Table 2:", nrow(tab2), "rows -> Table2_Main_Text.csv\n")

cat("\n== Self-check: Main-text Table 2 key genes ==\n")
chk <- c("GSTO2", "RNASEH2B", "FENDRR", "MTHFSD", "PPP1R14A",
         "CHMP2B", "MDM4", "IRX4", "SIDT1", "NKX3-1", "BIK", "CCND1")
print(tab2[Gene %in% chk])

writeLines(CAPTION_TEXT, file.path(OUT_DIR, "Table2_Caption.txt"))
cat("\n Caption written -> Table2_Caption.txt\n")

sessionInfo()
