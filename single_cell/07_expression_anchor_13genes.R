rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(dplyr)

if (!dir.exists("plots")) dir.create("plots")

bph_order12 <- c("T Cells", "B Cells", "Myeloid Cells", "Mast Cells", "Endothelial",
                 "Basal", "Club", "Luminal", "Fibroblasts", "Smooth Muscle")
pca_order12 <- c("T Cells", "B Cells", "Plasma cells", "Myeloid Cells", "Mast Cells",
                 "Endothelial", "Basal", "Club", "Luminal", "Malignant",
                 "Smooth Muscle", "Fibroblasts")

match_gene <- function(obj, g) {
  rn <- rownames(obj)
  hit <- grep(paste0("^", g, "$"), rn, ignore.case = TRUE, value = TRUE)
  if (length(hit) == 1) return(hit)
  hit2 <- grep(g, rn, ignore.case = TRUE, value = TRUE)
  if (length(hit2) >= 1) return(hit2[order(nchar(hit2))][1])
  return(NULL)
}

detect_groups <- function(obj, disease_keyword = "BPH|Tumor|Cancer") {
  grps <- sort(unique(as.character(obj$group)))
  is_dis <- grepl(disease_keyword, grps, ignore.case = TRUE)
  c(grps[!is_dis], grps[is_dis])
}

cat("Reading data...\n")
scRNA_bph <- readRDS("bph/newBPH.rds")
scRNA_pca <- readRDS("pca/newPCa.rds")

bph_map <- c("T Cells" = "T Cells", "B Cells" = "B Cells",
             "Myeloid Cells" = "Myeloid Cells", "Mast Cells" = "Mast Cells",
             "Endothelial" = "Endothelial", "Basal" = "Basal",
             "Club" = "Club", "Hillock" = "Club",
             "Luminal_AR+" = "Luminal", "Luminal_AR-" = "Luminal",
             "Pericytes" = "Smooth Muscle",
             "interstitial fibroblast" = "Fibroblasts",
             "peri-epithelial fibroblast" = "Fibroblasts",
             "vascular SMC(vSMC)" = "Smooth Muscle",
             "prostate SMC(pSMC)" = "Smooth Muscle")

pca_map <- c("T Cells" = "T Cells", "B Cells" = "B Cells",
             "Plasma cells" = "Plasma cells",
             "Myeloid Cells" = "Myeloid Cells", "Mast Cells" = "Mast Cells",
             "Endothelial" = "Endothelial", "Basal" = "Basal",
             "Luminal" = "Luminal", "Club" = "Club", "Malignant" = "Malignant",
             "SMC" = "Smooth Muscle", "Pericyte" = "Smooth Muscle",
             "iCAF" = "Fibroblasts", "myCAF" = "Fibroblasts", "apCAF" = "Fibroblasts")

scRNA_bph$cell_type_12ct <- factor(unname(bph_map[as.character(scRNA_bph$cell_type_sub)]),
                                   levels = bph_order12)
scRNA_pca$cell_type_12ct <- factor(unname(pca_map[as.character(scRNA_pca$cell_type_sub)]),
                                   levels = pca_order12)

bph_grps <- detect_groups(scRNA_bph, "BPH")
pca_grps <- detect_groups(scRNA_pca, "Tumor")
cat("BPH groups:", paste(bph_grps, collapse = ", "), "\n")
cat("PCa groups:", paste(pca_grps, collapse = ", "), "\n")

scRNA_bph <- subset(scRNA_bph, group %in% bph_grps[grepl("BPH", bph_grps, ignore.case = TRUE)])
scRNA_pca <- subset(scRNA_pca, group %in% pca_grps[grepl("Tumor", pca_grps, ignore.case = TRUE)])

bph_genes <- c("FENDRR", "MTHFSD", "RNASEH2B")
pca_genes <- c("CHMP2B", "PPP1R14A", "SIDT1", "C19ORF33", "MDM4",
               "BHLHA15", "SPINT2", "SLC22A2", "TPCN2", "WDPCP")

make_long <- function(obj, genes) {
  gene_real <- character()
  for (g in genes) {
    gr <- match_gene(obj, g)
    if (is.null(gr)) { cat("[ skipping ]", g, "Not found \n"); next }
    if (gr != g) cat("genes Info:", g, "->", gr, "\n")
    gene_real <- c(gene_real, gr)
    names(gene_real)[length(gene_real)] <- g
  }
  if (length(gene_real) == 0) return(NULL)

  expr_mat <- FetchData(obj, vars = gene_real, layer = "data")
  ct_vec <- as.character(obj$cell_type_12ct)

  long_df <- do.call(rbind, lapply(names(gene_real), function(g) {
    data.frame(ct = ct_vec, gene = g, expr = expr_mat[, gene_real[g]],
               stringsAsFactors = FALSE)
  }))
  long_df <- long_df[!is.na(long_df$ct), ]
  return(long_df)
}

cat("Info BPH Info...\n")
long_bph <- make_long(scRNA_bph, bph_genes)
cat("Info PCa Info...\n")
long_pca <- make_long(scRNA_pca, pca_genes)

keep_bph <- intersect(bph_genes, unique(long_bph$gene))
keep_pca <- intersect(pca_genes, unique(long_pca$gene))

make_tables <- function(long_df, genes_sub, order_ct, dataset_name, n_cells_vec) {
  if (is.null(long_df) || length(genes_sub) == 0) stop(dataset_name, "None Info genes")

  stat_raw <- long_df %>%
    dplyr::filter(gene %in% genes_sub) %>%
    dplyr::group_by(gene, ct) %>%
    dplyr::summarise(mean_expr = mean(expr),
                     pct_expr  = 100 * mean(expr > 0),
                     .groups = "drop")

  pad <- expand.grid(gene = genes_sub, ct = order_ct, stringsAsFactors = FALSE)
  stat <- pad %>%
    dplyr::left_join(stat_raw, by = c("gene", "ct")) %>%
    dplyr::mutate(gene = factor(gene, levels = genes_sub),
                  ct   = factor(ct, levels = order_ct)) %>%
    dplyr::arrange(gene, ct) %>%
    dplyr::mutate(dataset   = dataset_name,
                  n_cells   = as.integer(n_cells_vec[as.character(ct)]),
                  mean_expr = round(ifelse(is.na(mean_expr), 0, mean_expr), 3),
                  pct_expr  = round(ifelse(is.na(pct_expr), 0, pct_expr), 1))

  stat <- as.data.frame(stat)[, c("dataset", "gene", "ct",
                                  "n_cells", "mean_expr", "pct_expr")]
  names(stat)[names(stat) == "ct"] <- "cell_type"
  stat$gene <- as.character(stat$gene)
  stat$cell_type <- as.character(stat$cell_type)
  rownames(stat) <- NULL

  top_all <- stat %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  top_all <- as.data.frame(top_all)
  top_all$anchor_status <- ifelse(top_all$pct_expr >= 10, "anchored",
                                  ifelse(top_all$pct_expr > 0, "below_floor", "not_expressed"))
  top_all <- top_all[order(match(top_all$gene, genes_sub)), ]
  rownames(top_all) <- NULL

  cat("[", dataset_name, "] Info / Info:\n", sep = "")
  for (i in seq_len(nrow(top_all))) {
    cat("  ", top_all$gene[i], "->", top_all$cell_type[i],
        sprintf("(mean=%.2f, pct=%.0f%%, %s)\n",
                top_all$mean_expr[i], top_all$pct_expr[i], top_all$anchor_status[i]))
  }

  list(stat = stat, anchor = top_all)
}

n_bph <- table(scRNA_bph$cell_type_12ct)
n_pca <- table(scRNA_pca$cell_type_12ct)

tb_bph <- make_tables(long_bph, keep_bph, bph_order12, "BPH", n_bph)
tb_pca <- make_tables(long_pca, keep_pca, pca_order12, "PCa", n_pca)

tab_long   <- rbind(tb_bph$stat, tb_pca$stat)
tab_anchor <- rbind(tb_bph$anchor, tb_pca$anchor)

f_long   <- file.path("plots", "Table_S_gene_expr_12ct.csv")
f_anchor <- file.path("plots", "Table_S_gene_anchor.csv")
write.csv(tab_long,   f_long,   row.names = FALSE)
write.csv(tab_anchor, f_anchor, row.names = FALSE)

if (requireNamespace("writexl", quietly = TRUE)) {
  f_xlsx <- file.path("plots", "Table_S_gene_expression_12ct.xlsx")
  writexl::write_xlsx(list(Expr_full = tab_long, Anchor = tab_anchor), f_xlsx)
  cat("xlsx Written:", f_xlsx, "\n")
}

cat("\n=== Info Self-check ===\n")
for (fl in c(5, 10, 15, 20)) {
  a_bph <- tb_bph$stat %>%
    dplyr::filter(pct_expr >= fl) %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  a_pca <- tb_pca$stat %>%
    dplyr::filter(pct_expr >= fl) %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  cat(sprintf("floor=%d%%: BPH %d anchored, PCa %d anchored\n",
              fl, nrow(a_bph), nrow(a_pca)))
  print(as.data.frame(a_bph[, c("gene", "cell_type")]))
  print(as.data.frame(a_pca[, c("gene", "cell_type")]))
}

cat("Completed:\n", f_long, "\n  ", f_anchor,
    "\n(Full matrix", nrow(tab_long), "rows; Info", nrow(tab_anchor),
    "rows; disease group only)\n", sep = "")

sessionInfo()
