rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(166)

library(Seurat)

pca_epi_file <- "03_Epithelial_Annotation/epi_all_annotated.rds"
bph_epi_file <- "Epithelial_bph_subcluster_res0.6.rds"
is_pca <- file.exists(pca_epi_file)
is_bph <- file.exists(bph_epi_file)
if (is_pca == is_bph) {
  stop("Run this script from single_cell/pca or single_cell/bph with exactly one epithelial input available.")
}

group_name <- if (is_pca) "PCa" else "BPH"
epi_file <- if (is_pca) pca_epi_file else bph_epi_file
epi_col <- if (is_pca) "epi_subtype" else "cell_type"
output_file <- if (is_pca) "newPCa.rds" else "newBPH.rds"
big_file <- "scRNA_harmony_annotated.rds"
fib_files <- list.files(".", pattern = "^Fibroblast_.*final[.]rds$")
if (!file.exists(big_file)) stop("Missing input: ", big_file)
if (length(fib_files) != 1L) {
  print(fib_files)
  stop("Expected exactly one Fibroblast_*_final.rds file in the working directory.")
}
fib_file <- fib_files[1]

scRNA <- readRDS(big_file)
epi <- readRDS(epi_file)
fib <- readRDS(fib_file)
objects <- list(scRNA = scRNA, epi = epi, fib = fib)
for (object_name in names(objects)) {
  obj <- objects[[object_name]]
  cat(object_name, ":", nrow(obj), "genes,", ncol(obj), "cells\n")
  if (anyDuplicated(colnames(obj))) stop("Duplicated cell barcodes in ", object_name)
  if (anyDuplicated(colnames(obj@meta.data))) stop("Duplicated metadata columns in ", object_name)
}
if (!"cell_type" %in% colnames(scRNA@meta.data)) stop("Missing cell_type in the large object.")
if (!epi_col %in% colnames(epi@meta.data)) stop("Missing ", epi_col, " in the epithelial object.")
if (!"subcluster" %in% colnames(fib@meta.data)) stop("Missing subcluster in the fibroblast object.")
cat("Large-object cell types:\n")
print(table(scRNA$cell_type, useNA = "ifany"))
cat("Epithelial subtypes:\n")
print(table(epi@meta.data[[epi_col]], useNA = "ifany"))
cat("Fibroblast subtypes:\n")
print(table(fib$subcluster, useNA = "ifany"))

big_cells <- colnames(scRNA)
epi_hit <- intersect(colnames(epi), big_cells)
fib_hit <- intersect(colnames(fib), big_cells)
epi_rate <- length(epi_hit) / ncol(epi)
fib_rate <- length(fib_hit) / ncol(fib)
cat("Epithelial barcode match:", length(epi_hit), "/", ncol(epi), sprintf("(%.1f%%)\n", 100 * epi_rate))
cat("Fibroblast barcode match:", length(fib_hit), "/", ncol(fib), sprintf("(%.1f%%)\n", 100 * fib_rate))
if (epi_rate < 0.5 || fib_rate < 0.5) stop("Barcode match rate is below 50%; check barcode renaming.")
if (length(intersect(epi_hit, fib_hit)) > 0L) stop("Epithelial and fibroblast objects share cell barcodes.")

scRNA$cell_type_sub <- as.character(scRNA$cell_type)
epi_big_pos <- match(epi_hit, big_cells)
epi_sub_pos <- match(epi_hit, colnames(epi))
fib_big_pos <- match(fib_hit, big_cells)
fib_sub_pos <- match(fib_hit, colnames(fib))
scRNA$cell_type_sub[epi_big_pos] <- as.character(epi@meta.data[[epi_col]])[epi_sub_pos]
scRNA$cell_type_sub[fib_big_pos] <- as.character(fib$subcluster)[fib_sub_pos]
if (any(is.na(scRNA$cell_type_sub) | scRNA$cell_type_sub == "")) {
  stop("Missing subtype labels after barcode matching.")
}

parent_types <- c("Epithelial", "Fibroblasts", "Smooth Muscle")
drop_pos <- which(scRNA$cell_type_sub %in% parent_types)
cat("Cells not covered by subcluster annotation:", length(drop_pos), "\n")
if (length(drop_pos) > 0L) {
  print(table(scRNA$cell_type[drop_pos], useNA = "ifany"))
  scRNA <- subset(scRNA, cells = colnames(scRNA)[-drop_pos])
}
cat("Merged object:", nrow(scRNA), "genes,", ncol(scRNA), "cells\n")
print(table(scRNA$cell_type_sub, useNA = "ifany"))

if (!"harmony" %in% names(scRNA@reductions)) stop("The harmony reduction is missing.")
dims_use <- 1:10
if (max(dims_use) > ncol(Embeddings(scRNA, reduction = "harmony"))) {
  stop("The harmony reduction has fewer than 10 dimensions.")
}
scRNA <- RunUMAP(
  scRNA,
  reduction = "harmony",
  dims = dims_use,
  reduction.name = "umap_new",
  reduction.key = "umap_new_",
  verbose = FALSE
)
Idents(scRNA) <- "cell_type_sub"
saveRDS(scRNA, output_file)
stats <- as.data.frame(table(scRNA$cell_type_sub), stringsAsFactors = FALSE)
colnames(stats) <- c("cell_type_sub", "n_cells")
write.csv(stats, "cell_type_sub_stats.csv", row.names = FALSE)
cat("Saved:", output_file, "and cell_type_sub_stats.csv\n")
sessionInfo()
