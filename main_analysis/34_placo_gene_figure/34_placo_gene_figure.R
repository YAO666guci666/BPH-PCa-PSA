rm(list = ls())
options(stringsAsFactors = FALSE)
options(timeout = 120)

while (sink.number() > 0) sink()
while (dev.cur() > 1) dev.off()

if (!dir.exists("./Figure5")) dir.create("./Figure5", recursive = TRUE)
if (!dir.exists("./Figure5"))
  stop("Cannot create./Figure5, check write permissions:", getwd())

log_file <- "./Figure5/run_log.txt"
if (file.exists(log_file)) try(unlink(log_file, force = TRUE), silent = TRUE)
log_ok <- tryCatch({ sink(log_file, split = TRUE); TRUE },
                   error = function(e) {
                     cat("!! Cannot open the log file (", conditionMessage(e),
                         "), using console output only \n")
                     FALSE
                   })

cat("===== VER: Fig5_LAVA_v9_32_20261103 =====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("Working directory:", getwd(), "\n")
cat("Directory files:\n"); print(list.files("."))

library(data.table)
library(ggplot2)
library(ggrepel)
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n")

if (!requireNamespace("patchwork", quietly = TRUE)) {
  cat("First run: Installing patchwork(mirror)...\n")
  install.packages("patchwork",
                   repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
suppressPackageStartupMessages(library(patchwork))

col_bp   <- "#7C9895"
col_cp   <- "#1B3E76"
col_star <- "#FED71A"
pairs    <- c("BPH_PSA", "PCA_PSA")
pair_lab <- c("BPH–PSA", "PCa–PSA")
pair_col <- c("BPH–PSA" = col_bp, "PCa–PSA" = col_cp)
CUR_VER  <- "Fig5_LAVA_v9_32"

theme_hess <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(text = element_text(family = "serif", color = "black"),
          panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.6),
          panel.grid     = element_blank(),
          axis.line      = element_blank(),
          axis.ticks     = element_line(color = "black", linewidth = 0.4),
          axis.text      = element_text(size = base_size - 2.5),
          axis.title     = element_text(size = base_size - 0.5),
          plot.subtitle  = element_text(hjust = 0.5, size = base_size - 2,
                                        color = "grey35",
                                        margin = ggplot2::margin(b = 4)),
          strip.background = element_blank(),
          strip.text     = element_text(size = base_size - 1.5, face = "plain"),
          legend.position = "top",
          legend.title   = element_blank(),
          legend.text    = element_text(size = base_size - 2.5),
          legend.key.size = grid::unit(0.45, "cm"),
          legend.background = element_blank(),
          plot.margin    = ggplot2::margin(6, 10, 6, 10))
}

save_plot <- function(p, fbase, w, h) {
  bv <- attr(p, "build_ver")
  if (is.null(bv) || bv != CUR_VER)
    stop("\n The plot object came from an old script or a stale object (Info =", if (is.null(bv)) "None" else bv, ")\n",
         "The figure will not be saved. Restart the R session and rerun the complete script: Session -> Restart R, Info source(\"34_placo_gene_figure.R\") Info Running")
  for (f in paste0(fbase, c(".pdf", ".png"))) try(unlink(f, force = TRUE), silent = TRUE)
  use_base <- fbase
  pdf_ok <- tryCatch({
    suppressWarnings(ggsave(paste0(use_base, ".pdf"), p, width = w, height = h))
    TRUE
  }, error = function(e) {
    cat("!! PDF write failed:", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  if (!pdf_ok) {
    use_base <- paste0(fbase, "_", format(Sys.time(), "%H%M%S"))
    cat("File is locked, saving instead:", use_base, ".pdf\n", sep = "")
    suppressWarnings(ggsave(paste0(use_base, ".pdf"), p, width = w, height = h))
  }
  suppressWarnings(ggsave(paste0(use_base, ".png"), p, width = w, height = h,
                          dpi = 300, bg = "white"))
  cat("Written:", use_base, ".pdf / .png\n", sep = "")
  invisible(use_base)
}

find_fuma_file <- function(p, fname) {
  f1 <- file.path("../31_placo/PLACO_pairwise", p, p, fname)
  f2 <- file.path("../31_placo/PLACO_pairwise", p, fname)
  if (file.exists(f1)) return(f1)
  if (file.exists(f2)) return(f2)
  return(NULL)
}

cat("\n==== [0/3] Preparing liftOver chain ====\n")

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  cat("First run: Installing rtracklayer(mirror)...\n")
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
          BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("rtracklayer", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages(library(rtracklayer))

chain_candidates <- c(file.path(path.expand("~"), "liftchain",
                                "hg19ToHg38.over.chain.gz"),
                      file.path("./Figure5", "hg19ToHg38.over.chain.gz"))
cat("chain file candidate paths:\n"); print(chain_candidates)
cat("Existence by path:", paste(file.exists(chain_candidates), collapse = ", "), "\n")
chain_f <- chain_candidates[file.exists(chain_candidates)][1]
if (is.na(chain_f))
  stop("Not found in either location chain file. WSL Run in the terminal:\n",
       "  wget -O ~/liftchain/hg19ToHg38.over.chain.gz \"https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz\"\n",
       "then rerun this script")
cat("Using chain file:", normalizePath(chain_f), " (", file.info(chain_f)$size, "bytes)\n", sep = "")

chain_ok <- function(f) {
  tryCatch({
    con <- gzcon(file(f, "rb"))
    on.exit(close(con))
    ln <- readLines(con, n = 1, warn = FALSE)
    length(strsplit(ln, "\\s+")[[1]]) >= 11
  }, error = function(e) FALSE)
}
if (!chain_ok(chain_f))
  stop("chain file corrupted (insufficient fields in the first line 11):", chain_f, "-- download it again before rerunning")
cat("chain file Content validation passed, Decompressing and importing...\n")
plain_f <- file.path(tempdir(), "hg19ToHg38.over.chain")
con_in  <- gzcon(file(chain_f, "rb"))
lines   <- readLines(con_in, warn = FALSE)
close(con_in)
writeLines(lines, plain_f)
chain <- import.chain(plain_f)
cat("chain Import successful \n")

can_gr   <- GRanges(seqnames = "chr6",
                    ranges = IRanges(start = 160581374, end = 160581374))
can_lift <- liftOver(can_gr, chain)
if (length(can_lift[[1]]) == 0) stop("liftOver Self-check failed: Info could not be mapped")
cr <- range(can_lift[[1]])
if (as.character(seqnames(cr)) != "chr6" || start(cr) > 160160342 || end(cr) < 160160342)
  stop("liftOver Self-check failed: Info expected mapping to 160160342, Info",
       as.character(seqnames(cr)), ":", start(cr), "-", end(cr), "(incorrect chain file?)")
cat("Info Self-check passed: 6:160581374(hg19) ->", as.character(seqnames(cr)), ":",
    start(cr), "-", end(cr), " (GRCh38)\n", sep = "")

cat("\n==== [1/3] lift leadSNPs hg19 -> GRCh38 (Info) ====\n")

for (p in pairs) {
  out_f <- file.path("./Figure5", paste0("leadSNPs_", p, ".GRCh38.txt"))
  aud_f <- file.path("./Figure5", paste0("leadSNPs_", p, ".audit.csv"))
  if (file.exists(out_f) && file.exists(aud_f)) {
    amed <- tryCatch(fread(aud_f)[, median(abs(shift_bp))], error = function(e) NA_real_)
    if (!is.na(amed) && amed > 1000) {
      cat("  [", p, "] Using the completed file (median shift", round(amed / 1000, 1),
          "kb, not repeated lift):", out_f, "\n", sep = "")
      next
    }
    cat("  [", p, "] Existing output failed validation, rerunning lift\n", sep = "")
  }
  src <- find_fuma_file(p, "leadSNPs.txt")
  if (is.null(src)) stop("Not found Source file:../31_placo/PLACO_pairwise/", p, "/leadSNPs.txt")
  cat("  [", p, "] Source (read-only):", src, "\n", sep = "")
  d <- fread(src)
  n <- nrow(d)
  d[, chr := suppressWarnings(as.integer(chr))]
  if (any(is.na(d$chr)) || any(!d$chr %in% 1:22))
    stop(p, "Contains non- 1-22 chromosomes rows, manual processing is required")

  gr     <- GRanges(seqnames = paste0("chr", d$chr),
                    ranges = IRanges(start = d$pos, end = d$pos))
  lifted <- liftOver(gr, chain)
  nmul   <- elementNROWS(lifted)

  chr_v <- rep(NA_integer_, n); pos_v <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    L <- lifted[[i]]
    if (length(L) == 0) next
    L <- L[as.character(seqnames(L)) != "*"]
    if (length(L) == 0) next
    j <- which.max(width(L))
    chr_v[i] <- as.integer(sub("^chr", "", as.character(seqnames(L))[j]))
    pos_v[i] <- start(L)[j]
  }
  if (any(nmul > 1))
    cat("  !! ", sum(nmul > 1), "SNP Info, Info (Info nmulti columns)\n", sep = "")

  bad <- which(is.na(chr_v) | chr_v != d$chr)
  if (length(bad) > 0) {
    fcf <- file.path("./Figure5", paste0("lead_failures_", p, ".csv"))
    fwrite(data.table(row = bad, rsID = d$rsID[bad], chr = d$chr[bad],
                      pos = d$pos[bad]), fcf, sep = "\t")
    stop("[", p, "] ", length(bad), "SNP None Info, Info:", fcf,
         "\n manual processing is required")
  }

  d[, chr_hg19 := chr][, pos_hg19 := pos]
  d[, `:=`(chr = chr_v, pos = pos_v)]
  ua <- strsplit(d$uniqID, ":")
  d[, uniqID := paste0(chr, ":", pos, ":",
                       vapply(ua, function(x) x[3], character(1)), ":",
                       vapply(ua, function(x) x[4], character(1)))]
  fwrite(d, out_f, sep = "\t")
  fwrite(data.table(rsID = d$rsID, chr_hg19 = d$chr_hg19, pos_hg19 = d$pos_hg19,
                    chr_38 = d$chr, pos_38 = d$pos,
                    shift_bp = d$pos - d$pos_hg19, nmulti = nmul), aud_f, sep = "\t")
  cat("  [", p, "] Completed", n, "rows; median shift",
      round(median(d$pos - d$pos_hg19) / 1000, 1),
      "kb\n Written:", out_f, " + ", aud_f, "\n", sep = "")
  if (p == "PCA_PSA") {
    chk <- fread(out_f)
    rr  <- chk[tolower(rsID) == "rs651164"]
    if (nrow(rr) != 1 || rr$pos != 160160342)
      stop("chain Info: rs651164 Info 160160342, Info",
           if (nrow(rr)) rr$pos else "Not found", "-- check the log")
  }
}

cat("\n==== [2/3] lift LAVA Info hg19 -> GRCh38 (Info) ====\n")

manual_fix <- data.table(pair = "BPH_PSA", row = 266, chr = 18,
                         start = 15220936, stop = 22429733)

lava_targets <- c(BPH_PSA = "../33_lava/PSA_BPH/PSA_BPH_all_result.txt", PCA_PSA = "../33_lava/PSA_PC/PSA_PC_all_result.txt")
for (nm in names(lava_targets)) {
  src <- lava_targets[nm]
  out_f <- file.path("./Figure5", sub("\\.txt$", ".GRCh38.txt", basename(src)))
  aud_f <- file.path("./Figure5", sub("all_result\\.txt$", "lift_audit.csv", basename(src)))
  if (file.exists(out_f) && file.exists(aud_f)) {
    amed <- tryCatch(fread(aud_f)[, median(abs(shift_kb))], error = function(e) NA_real_)
    if (!is.na(amed) && amed > 100) {
      cat("  [", nm, "] Using the completed file (median shift", round(amed, 1),
          "kb, not repeated lift):", out_f, "\n", sep = "")
      next
    }
    cat("  [", nm, "] Info Output Info median shift Info (", round(amed, 1), "kb), rerunning lift\n", sep = "")
  }
  if (!file.exists(src)) stop("Not found Source file:", src)
  cat("  [", nm, "] Source (read-only):", src, "\n", sep = "")
  d <- fread(src)
  setnames(d, tolower(names(d)))
  if (!all(c("chr", "start", "stop") %in% names(d)))
    stop(src, "Missing chr/start/stop columns")
  d[, chr := suppressWarnings(as.integer(chr))]
  if (any(is.na(d$chr)) || any(!d$chr %in% 1:22))
    stop(src, "Contains non- 1-22 chromosomes rows, manual processing is required")
  n <- nrow(d)

  gr     <- GRanges(seqnames = paste0("chr", d$chr),
                    ranges = IRanges(start = d$start, end = d$stop))
  lifted <- liftOver(gr, chain)
  pieces <- elementNROWS(lifted)

  chr_v <- rep(NA_integer_, n); st_v <- rep(NA_real_, n); en_v <- rep(NA_real_, n)
  cov_v <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    L <- lifted[[i]]
    if (length(L) == 0) next
    L <- L[as.character(seqnames(L)) != "*"]
    if (length(L) == 0) next
    cov  <- tapply(width(L), as.character(seqnames(L)), sum)
    best <- names(cov)[which.max(cov)]
    Lb   <- L[as.character(seqnames(L)) == best]
    r    <- range(Lb)
    chr_v[i] <- as.integer(sub("^chr", "", best))
    st_v[i]  <- start(r)
    en_v[i]  <- end(r)
    cov_v[i] <- sum(width(Lb)) / (d$stop[i] - d$start[i] + 1)
  }
  if (any(pieces > 1))
    cat("  !! ", sum(pieces > 1), "Info, Info chromosomes Merging \n", sep = "")

  mf <- manual_fix[pair == nm]
  if (nrow(mf) > 0) {
    for (k in seq_len(nrow(mf))) {
      m <- mf[k]
      if (m$row < 1 || m$row > n) stop("manual_fix rows Info:", m$row, " (", nm, ")")
      chr_v[m$row] <- m$chr; st_v[m$row] <- m$start; en_v[m$row] <- m$stop
      cov_v[m$row] <- 1
    }
    cat("Info", nrow(mf), "Info (", nm, ")\n", sep = "")
  }

  bad <- which(is.na(chr_v) | cov_v < 0.5 | chr_v != d$chr)
  if (length(bad) > 0) {
    fcf <- file.path("./Figure5", paste0("lift_failures_", nm, ".csv"))
    fwrite(data.table(row = bad, chr = d$chr[bad], start = d$start[bad],
                      stop = d$stop[bad], cov = round(cov_v[bad], 2)), fcf, sep = "\t")
    stop("[", nm, "] ", length(bad), "Info None Info, Info:", fcf,
         "\n Inspect this file, Add the corrected coordinates to manual_fix and rerun")
  }

  d[, chr_hg19 := chr][, start_hg19 := start][, stop_hg19 := stop]
  d[, `:=`(chr = chr_v, start = st_v, stop = en_v)]
  fwrite(d, out_f, sep = "\t")
  fwrite(data.table(row = seq_len(n), chr_hg19 = d$chr_hg19, start_hg19 = d$start_hg19,
                    stop_hg19 = d$stop_hg19, chr_38 = d$chr,
                    start_38 = d$start, stop_38 = d$stop,
                    shift_kb = round((d$start - d$start_hg19) / 1000, 1),
                    pieces = pieces, cov = round(cov_v, 3)), aud_f, sep = "\t")
  cat("  [", nm, "] Completed", n, "Info; median shift",
      round(median(d$start - d$start_hg19) / 1000, 1),
      "kb\n Written:", out_f, " + ", aud_f, "\n", sep = "")
}

cat("\n==== [3/3] Info ====\n")

read_lead <- function(p) {
  f <- file.path("./Figure5", paste0("leadSNPs_", p, ".GRCh38.txt"))
  if (!file.exists(f))
    stop("Not found", f, "(1 Info Generating; Info, Info Running Info Script)")
  cat("  [", p, "] Reading:", f, " (", nrow(fread(f)), "rows)\n", sep = "")
  d <- fread(f); setnames(d, tolower(names(d)))
  if (!"rsid" %in% names(d)) setnames(d, names(d)[1], "rsid")
  if ("p" %in% names(d) && !"p_placo" %in% names(d)) setnames(d, "p", "p_placo")
  d[, `:=`(pair = p, pair_key = p)]
  d[, .(rsid, chr = as.integer(chr), pos = as.numeric(pos),
        p_placo = as.numeric(p_placo), pair, pair_key)]
}

lava_files <- c(BPH_PSA = "./Figure5/PSA_BPH_all_result.GRCh38.txt",
                PCA_PSA = "./Figure5/PSA_PC_all_result.GRCh38.txt")
lava <- rbindlist(lapply(names(lava_files), function(nm){
  f <- lava_files[nm]
  if (!file.exists(f)) stop("Not found", f, "(2 Info Generating; Info, Info Running Info Script)")
  d <- fread(f)
  setnames(d, tolower(names(d)))
  need <- c("chr", "start", "stop", "rho", "p")
  miss <- setdiff(need, names(d))
  if (length(miss) > 0)
    stop(nm, "Missing columns:", paste(miss, collapse = ", "),
         "; Actual columns:", paste(names(d), collapse = ", "))
  d[, .(pair_key = nm, chr = as.integer(chr), start = as.numeric(start),
        stop = as.numeric(stop), rho = as.numeric(rho), p = as.numeric(p))]
}))
lava <- lava[!is.na(p)]
lava[, ntest := .N, by = pair_key]
lava[, fdr  := p.adjust(p, method = "BH"), by = pair_key]
lava[, sig  := fdr < 0.05]
lava[, logp := -log10(p)]
cat("LAVA Info:",
    lava[, .N, by = pair_key][, paste(pair_key, N, sep = "=", collapse = ", ")],
    "Significant (FDR<0.05):",
    lava[, sum(sig), by = pair_key][, paste(pair_key, V1, sep = "=", collapse = ", ")], "\n")

lead <- rbindlist(lapply(pairs, read_lead))
lead[, pair := factor(pair, levels = pairs, labels = pair_lab)]
lead[, logp := -log10(p_placo)]
lead <- lead[chr %in% 1:22]

ext <- rbind(lava[, .(chr, s = start, e = stop)],
             lead[, .(chr, s = pos, e = pos)])
chr_range <- ext[, .(xmin = min(s), xmax = max(e)), by = chr]

map_band <- function(chr_v, pos_v) {
  x   <- as.numeric(chr_v)
  xmn <- chr_range$xmin[match(chr_v, chr_range$chr)]
  xmx <- chr_range$xmax[match(chr_v, chr_range$chr)]
  ok  <- !is.na(xmn) & (xmx > xmn)
  x[ok] <- chr_v[ok] - 0.36 + 0.72 * (pos_v[ok] - xmn[ok]) / (xmx[ok] - xmn[ok])
  x
}

lava[, `:=`(x_s = map_band(chr, start), x_e = map_band(chr, stop))]
lead[, x := map_band(chr, pos)]
lead[, lid := .I]

lead[, in_lava := FALSE]
for (i in seq_len(nrow(lava))) {
  if (!lava$sig[i]) next
  lead[pair_key == lava$pair_key[i] & chr == lava$chr[i] &
         pos >= lava$start[i] & pos <= lava$stop[i], in_lava := TRUE]
}
cat("LAVA Significant Info lead SNP:", lead[, sum(in_lava)], "\n")
if (lead[, sum(in_lava)] == 0)
  stop("in_lava Info FALSE: LAVA and lead file build Info? Check./Figure5 GRCh38 file")

fig_genes <- c("C19ORF33","CHMP2B","MDM4","PPP1R14A","BHLHA15","SIDT1","SLC22A2",
               "SPINT2","TPCN2","WDPCP","FENDRR","MTHFSD","RNASEH2B")
gene2pair <- c(C19ORF33 = "PCA_PSA", CHMP2B = "PCA_PSA", MDM4 = "PCA_PSA",
               PPP1R14A = "PCA_PSA", BHLHA15 = "PCA_PSA", SIDT1 = "PCA_PSA",
               SLC22A2 = "PCA_PSA", SPINT2 = "PCA_PSA", TPCN2 = "PCA_PSA",
               WDPCP = "PCA_PSA", FENDRR = "BPH_PSA", MTHFSD = "BPH_PSA",
               RNASEH2B = "BPH_PSA")

rest38_one <- function(sym) {
  u <- paste0("https://rest.ensembl.org/lookup/symbol/homo_sapiens/",
              sym, "?content-type=application/json")
  for (wait in c(1, 5, 15, 30)) {
    txt <- tryCatch(paste(readLines(u, warn = FALSE), collapse = ""),
                    error = function(e) NULL)
    if (!is.null(txt) && grepl('"start":', txt)) {
      st <- sub('"start":', '', regmatches(txt, regexpr('"start":[0-9]+', txt)))
      en <- sub('"end":',   '', regmatches(txt, regexpr('"end":[0-9]+', txt)))
      sr <- gsub('"', '', regmatches(txt, regexpr('"seq_region_name":"[^"]+"', txt)))
      sr <- sub("seq_region_name:", '', sr)
      id <- gsub('"', '', regmatches(txt, regexpr('"id":"ENSG[0-9]+"', txt)))
      id <- sub("id:", '', id)
      return(list(chr = suppressWarnings(as.integer(sr)), start = as.numeric(st),
                  end = as.numeric(en), ens = id, sym_used = sym))
    }
    Sys.sleep(wait)
  }
  NULL
}
rest38 <- function(g) {
  cand <- unique(c(g, paste0(toupper(substr(g, 1, 1)), tolower(substring(g, 2)))))
  for (s in cand) {
    r <- rest38_one(s)
    if (!is.null(r)) return(r)
  }
  NULL
}

GC_FILE <- "./Figure5/gene_coord_GRCh38.tsv"
gene_coord <- NULL
cache_ok <- FALSE
if (file.exists(GC_FILE)) {
  gchk <- fread(GC_FILE)
  if (nrow(gchk[gene == "SLC22A2" & gstart > 160150000 & gstart < 160300000]) == 1 &&
      all(fig_genes %in% gchk$gene) &&
      gchk[gene %in% fig_genes, all(!is.na(gchr) & !is.na(gstart) & !is.na(gstop))]) {
    gene_coord <- gchk[gene %in% fig_genes]
    cache_ok <- TRUE
    cat("Using genes Info cache (Info", length(fig_genes), "genes, skipping REST):",
        GC_FILE, "\n", sep = "")
  } else {
    cat("genes Info cache Self-check Info, Info REST Info \n")
  }
}

if (!cache_ok) {
  cat("Info genes Info (Ensembl GRCh38 REST)...\n")
  can <- rest38_one("SLC22A2")
  if (is.null(can))
    stop("REST Info (Info): Check / proxy, or Info Running. \n",
         "After one successful run, cache", GC_FILE, "will be reused.")
  if (can$chr != 6 || can$start < 160150000 || can$start > 160300000)
    stop("Self-check failed: SLC22A2 GRCh38 Info 160.15-160.30 Mb, Info",
         can$start, "(REST unexpected build returned, do not use)")
  cat("Info Self-check passed: SLC22A2 = chr", can$chr, ":", can$start, "-", can$end, "\n", sep = "")

  gc_list <- lapply(fig_genes, rest38)
  names(gc_list) <- fig_genes
  bad <- fig_genes[sapply(gc_list, is.null)]
  if (length(bad) > 0)
    stop("GRCh38 REST Info:", paste(bad, collapse = ", "), "(Info Time Info Running)")
  fixed <- fig_genes[sapply(fig_genes, function(g) gc_list[[g]]$sym_used != g)]
  if (length(fixed) > 0)
    cat("Info:", paste(paste0(fixed, " -> ", sapply(fixed, function(g) gc_list[[g]]$sym_used)),
                           collapse = ", "), "\n", sep = "")
  new_coord <- rbindlist(lapply(fig_genes, function(g) {
    p <- gc_list[[g]]
    data.table(gene = g, ENSG_used = p$ens, gchr = p$chr,
               gstart = p$start, gstop = p$end, coord_source = "Ensembl GRCh38 REST")
  }))
  if (file.exists(GC_FILE)) {
    old_keep <- fread(GC_FILE)[!gene %in% new_coord$gene]
    gene_coord <- rbind(old_keep, new_coord, fill = TRUE)
  } else {
    gene_coord <- new_coord
  }
  fwrite(gene_coord[order(gene)], GC_FILE, sep = "\t")
}
cat("genes Info (Info Source):\n")
print(as.data.frame(gene_coord[gene %in% fig_genes][order(gene)]))

star_pairs <- rbindlist(lapply(fig_genes, function(g) {
  gc <- gene_coord[gene == g]
  pk <- gene2pair[[g]]
  lp <- lead[pair_key == pk & chr == gc$gchr]
  if (nrow(lp) == 0) return(NULL)
  d <- ifelse(lp$pos >= gc$gstart & lp$pos <= gc$gstop, 0,
              pmin(abs(lp$pos - gc$gstart), abs(lp$pos - gc$gstop)))
  j <- which.min(d)
  data.table(idx = lp$lid[j], gene = g, dist_kb = round(d[j] / 1000, 1))
}))
if (is.null(star_pairs)) star_pairs <- data.table(idx = integer(0), gene = character(0),
                                                  dist_kb = numeric(0))
star_pairs <- unique(star_pairs, by = "gene")
star_pairs[, pair_key := sapply(gene, function(g) gene2pair[[g]])]
cat("Info (kb, Info):\n")
print(as.data.frame(star_pairs[order(-dist_kb)]))
if (any(star_pairs$dist_kb > 1000))
  stop("Info >1Mb(Info LD Info), Info Output:",
       paste(star_pairs$gene[star_pairs$dist_kb > 1000], collapse = ", "))

anchor_audit <- merge(star_pairs,
                      gene_coord[, .(gene, ENSG_used, coord_source, gchr, gstart, gstop)],
                      by = "gene")
anchor_audit <- merge(anchor_audit,
                      lead[, .(lid, anchor_rsid = rsid, anchor_chr = chr, anchor_pos = pos)],
                      by.x = "idx", by.y = "lid")
anchor_audit[, pair := ifelse(pair_key == "BPH_PSA", "BPH–PSA", "PCa–PSA")]
fwrite(anchor_audit[order(gene)], "./Figure5/Table_gene_anchor_check.csv")
cat("Info:./Figure5/Table_gene_anchor_check.csv\n")

ii <- star_pairs$idx
sl <- data.table(
  pair     = lead$pair[ii],
  rsid     = as.character(lead$rsid[ii]),
  chr      = lead$chr[ii],
  pos      = lead$pos[ii],
  p_placo  = lead$p_placo[ii],
  logp     = lead$logp[ii],
  top_gene = as.character(star_pairs$gene)
)
star_lab <- sl[, .(genes = paste(sort(unique(top_gene)), collapse = " / "),
                   chr = chr[1], pos = pos[1],
                   logp = logp[which.min(p_placo)]),
               by = .(pair, rsid)]
star_lab[, lab := sprintf("%s · %s", genes, rsid)]
star_lab[, x := map_band(chr, pos)]
covered <- unique(unlist(strsplit(star_lab$genes, " / ")))
cat("Info:", length(covered), "genes /", nrow(star_lab), "Info (Info",
    length(fig_genes), "genes)\n", sep = "")
miss_g2 <- setdiff(fig_genes, covered)
if (length(miss_g2) > 0) cat("!! Info genes:", paste(miss_g2, collapse = ", "), "\n")
print(star_lab[, .(pair, genes, rsid, chr, pos)])

gene_regions <- rbindlist(lapply(fig_genes, function(g){
  gc <- gene_coord[gene == g]
  pk <- gene2pair[[g]]
  reg <- lava[sig == TRUE & pair_key == pk & chr == gc$gchr &
                start <= gc$gstop & stop >= gc$gstart]
  data.table(gene = g,
             pair = ifelse(pk == "BPH_PSA", "BPH–PSA", "PCa–PSA"),
             in_lava_region = nrow(reg) > 0,
             lava_rho = ifelse(nrow(reg) > 0, reg$rho[1], NA_real_),
             lava_p   = ifelse(nrow(reg) > 0, reg$p[1], NA_real_),
             lava_fdr = ifelse(nrow(reg) > 0, reg$fdr[1], NA_real_))
}))
cat("genes Info LAVA Significant Info: matched", sum(gene_regions$in_lava_region),
    "/", nrow(gene_regions), "\n")
print(as.data.frame(gene_regions))
fwrite(gene_regions, "./Figure5/Table_gene_LAVA_regions.csv")

lead[, is_toplocus := lid %in% star_pairs$idx]
stars <- lead[is_toplocus == TRUE]
hits  <- lead[in_lava == TRUE]

cat("  lead$pair:\n");     print(table(lead$pair, useNA = "ifany"))
cat("  stars$pair:\n");    print(table(stars$pair, useNA = "ifany"))
cat("  hits$pair:\n");     print(table(hits$pair, useNA = "ifany"))
if (anyNA(lead$pair) || anyNA(stars$pair) || anyNA(hits$pair))
  stop("pair Info NA")
cat("validating Info: lead", nrow(lead), "/ hits", nrow(hits),
    "/ stars", nrow(stars), "/ labelled loci", nrow(star_lab), "\n")

plab_of <- function(pk0) ifelse(pk0 == "BPH_PSA", "BPH–PSA", "PCa–PSA")

cat("\n---- Info 5A: combined ----\n")

dt_band <- lava[sig == TRUE, .(x_s, x_e, pair = plab_of(pair_key),
                               panel = factor(plab_of(pair_key), levels = pair_lab))]

dt_band_lab <- rbindlist(lapply(pairs, function(pk0) {
  d <- lava[pair_key == pk0 & sig == TRUE][order(p)]
  if (nrow(d) == 0) return(NULL)
  d <- d[seq_len(min(5, nrow(d)))]
  d[, .(x = (x_s + x_e) / 2,
        lab = sprintf("chr%d:%.2f\u2013%.2f Mb", chr, start / 1e6, stop / 1e6),
        panel = factor(plab_of(pk0), levels = pair_lab))]
}))
dt_band_lab <- dt_band_lab[order(panel, x)]
dt_band_lab[, w := 0.12 * nchar(lab) + 0.6]
lvl_v <- integer(nrow(dt_band_lab))
for (pp in unique(dt_band_lab$panel)) {
  idx <- which(dt_band_lab$panel == pp)
  lx <- -Inf; lw <- 0
  for (k in idx) {
    placed <- FALSE
    for (L in seq_along(lx)) {
      if (dt_band_lab$x[k] - lx[L] >= (dt_band_lab$w[k] + lw[L]) / 2 + 0.35) {
        lvl_v[k] <- L - 1L
        lx[L] <- dt_band_lab$x[k]; lw[L] <- dt_band_lab$w[k]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      lvl_v[k] <- length(lx)
      lx <- c(lx, dt_band_lab$x[k]); lw <- c(lw, dt_band_lab$w[k])
    }
  }
}
dt_band_lab[, lvl := lvl_v]
cat("Info (Info / Info /x/ Info / Info):\n")
print(as.data.frame(dt_band_lab[, .(panel, lab, x = round(x, 2), w = round(w, 2), lvl)]))
if (max(dt_band_lab$lvl) > 3)
  stop("Info 4 Info, Info Output: Info 2.6 Info 2.4 Info Running")

dt_dot_all <- lead[, .(x, y = logp, panel = pair)]
dt_dot_hit <- lead[in_lava == TRUE, .(x, y = logp, pair = plab_of(pair_key),
                                      panel = pair)]
dt_star <- stars[, .(x, y = logp, panel = pair)]
dt_lab  <- star_lab[, .(pair = as.character(pair), x, y = logp, lab)]
dt_lab[, panel := factor(pair, levels = pair_lab)]

ymax_all <- max(lead$logp, na.rm = TRUE)
dt_band_lab[, y := ymax_all * (1.10 + 0.065 * lvl)]

dt_thr <- data.table(panel = factor(pair_lab, levels = pair_lab),
                     ythr = -log10(5e-8), lab = "P = 5e-8")

pA <- ggplot(dt_dot_all, aes(x = x, y = y)) +
  geom_rect(data = dt_band, inherit.aes = FALSE,
            aes(xmin = x_s, xmax = x_e, fill = pair),
            ymin = -Inf, ymax = Inf, alpha = 0.20, color = NA) +
  geom_segment(data = dt_thr, inherit.aes = FALSE,
               aes(x = 0.62, xend = 22.42, y = ythr, yend = ythr),
               color = "black", linetype = "dotted", linewidth = 0.4) +
  geom_text(data = dt_thr, inherit.aes = FALSE,
            aes(x = 22.42, y = ythr, label = lab, hjust = 1, vjust = -0.5),
            family = "serif", size = 2.6, color = "grey40") +
  geom_point(data = dt_dot_all, inherit.aes = FALSE,
             aes(x = x, y = y, shape = "All lead SNPs"),
             color = "grey80", size = 1.4, alpha = 0.85) +
  geom_point(data = dt_dot_hit, inherit.aes = FALSE,
             aes(x = x, y = y, color = pair),
             shape = 16, size = 1.9, alpha = 0.95) +
  geom_point(data = dt_star, inherit.aes = FALSE,
             aes(x = x, y = y, shape = "High-score gene locus"),
             color = "black", fill = col_star, size = 3.4, stroke = 1.0) +
  geom_label(data = dt_band_lab, inherit.aes = FALSE,
             aes(x = x, y = y, label = lab),
             family = "serif", size = 2.6, fontface = "plain", color = "black",
             fill = alpha("white", 0.9), label.size = 0.25,
             label.padding = unit(0.15, "lines"), show.legend = FALSE) +
  geom_label_repel(data = dt_lab, inherit.aes = FALSE,
                   aes(x = x, y = y, label = lab),
                   ylim = c(0, ymax_all * 1.06),
                   family = "serif", size = 2.9, fontface = "plain", color = "black",
                   fill = alpha("white", 0.88), label.size = 0.25, label.padding = 0.22,
                   max.overlaps = Inf, box.padding = 0.55, point.padding = 0.4,
                   min.segment.length = 0.1, segment.color = "grey45",
                   segment.size = 0.3, seed = 42, show.legend = FALSE) +
  facet_wrap(~ panel, ncol = 2) +
  scale_x_continuous(breaks = 1:22, labels = 1:22,
                     limits = c(0.5, 22.5), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, ymax_all * 1.36), expand = c(0, 0)) +
  scale_fill_manual(values = pair_col, breaks = names(pair_col)) +
  scale_color_manual(values = pair_col, breaks = names(pair_col), guide = "none") +
  scale_shape_manual(values = c("All lead SNPs" = 16,
                                "High-score gene locus" = 23)) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
         shape = guide_legend(order = 2)) +
  labs(x = "Chromosome", y = "−log10 P-value (PLACO)") +
  theme_hess(12)
attr(pA, "build_ver") <- CUR_VER

msgA <- tryCatch({ ggplot_build(pA); NULL }, error = function(e) conditionMessage(e))
if (!is.null(msgA)) { print(msgA); stop("Info 5A Building Error, Info") }

cat("Info A(combined):\n")
save_plot(pA, "./Figure5/Figure5A_Shared_Loci", 12, 5.6)
cat("\n")

cat("\n---- Z-score concordance(Info -only) ----\n")

gwas_files <- list(BPH = "../01_preprocess_gwas/BPH.txt", PSA = "../01_preprocess_gwas/PSA.txt", PCA = "../01_preprocess_gwas/PCA.txt")
gwas <- lapply(names(gwas_files), function(nm){
  d <- fread(gwas_files[[nm]])
  setnames(d, tolower(names(d)))
  if ("snp" %in% names(d)) setnames(d, "snp", "SNP")
  if (!"SNP" %in% names(d) && "rsid" %in% names(d)) setnames(d, "rsid", "SNP")
  if ("beta" %in% names(d)) setnames(d, "beta", "BETA")
  if ("se" %in% names(d)) setnames(d, "se", "SE")
  if (!all(c("SNP", "BETA", "SE") %in% names(d)))
    stop(nm, "Missing SNP/rsID/BETA/SE columns, Actual columns:", paste(names(d), collapse = ", "))
  d[, .(SNP, Z = BETA / SE)]
})
names(gwas) <- names(gwas_files)

pair_traits <- list(`BPH–PSA` = c("BPH", "PSA"), `PCa–PSA` = c("PCA", "PSA"))

conc <- rbindlist(lapply(names(pair_traits), function(plab){
  tr <- pair_traits[[plab]]
  pdir <- ifelse(plab == "BPH–PSA", "BPH_PSA", "PCA_PSA")
  ld <- read_lead(pdir)
  m1 <- match(ld$rsid, gwas[[tr[1]]]$SNP)
  m2 <- match(ld$rsid, gwas[[tr[2]]]$SNP)
  cat(sprintf("[%s] Z Info: %s %d/%d, %s %d/%d\n", plab,
              tr[1], sum(!is.na(m1)), nrow(ld), tr[2], sum(!is.na(m2)), nrow(ld)))
  data.table(pair = plab, rsid = ld$rsid, Z1 = gwas[[tr[1]]]$Z[m1],
             Z2 = gwas[[tr[2]]]$Z[m2])
}))
conc <- conc[!is.na(Z1) & !is.na(Z2)]
conc[, same_dir := sign(Z1) == sign(Z2)]

stat <- conc[, .(n = .N,
                 concordant = sum(same_dir),
                 pct = round(100 * mean(same_dir), 1),
                 r = round(cor(Z1, Z2), 2),
                 slope = round(unname(coef(lm(Z1 ~ Z2))[2]), 2)), by = pair]
cat("=== Z concordance Info (Info Main-text Info, Info) ===\n"); print(stat)
fwrite(stat, "./Figure5/Table_Z_concordance_stats.csv")
cat("Written:./Figure5/Table_Z_concordance_stats.csv\n")
cat("Main-text Info: %d/%d lead SNPs showed concordant effect directions (r = %.2f)\n",
    stat$concordant[1], stat$n[1], stat$r[1])

cat("\n===== VER Fig5_LAVA_v9_32 Info rows Completed =====\n")
if (log_ok) sink()

sessionInfo()
