options(stringsAsFactors = FALSE)
rm(list=ls())
gc()

library(Seurat)
library(tidyverse)
library(dplyr)
library(harmony)
library(patchwork)

message("=== BPH Info cells Info (Info Harmony Info)===")

scRNA = read_rds('scRNA_harmony_annotated.rds')
Idents(scRNA) = scRNA@meta.data$cell_type

Epithelial_scRNA = subset(scRNA, subset = cell_type == "Epithelial")
rm(scRNA)
message(sprintf("Info cells: %d", ncol(Epithelial_scRNA)))

Epithelial_scRNA_1 = Epithelial_scRNA %>%
  NormalizeData()

message("\n Info...")
Epithelial_scRNA_1[["mt_percent"]] <- PercentageFeatureSet(Epithelial_scRNA_1, pattern = "^MT-")

before <- ncol(Epithelial_scRNA_1)
Epithelial_scRNA_1 <- subset(Epithelial_scRNA_1,
                             subset = nFeature_RNA > 200 &
                               nFeature_RNA < 8000 &
                               mt_percent < 20)
after <- ncol(Epithelial_scRNA_1)
message(sprintf("Info: %d > %d (Info %.1f%%)", before, after, 100*(before-after)/before))

Epithelial_scRNA_1 <- Epithelial_scRNA_1 %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

message("\n Generating Info PCA Info...")

p_before <- DimPlot(Epithelial_scRNA_1,
                    reduction = 'pca',
                    group.by = 'group',
                    pt.size = 0.5,
                    cols = c("BPH" = "#E64B35", "Healthy" = "#27AE60")) +
  ggtitle("Before Harmony: PCA by Group") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

p_before_sample <- DimPlot(Epithelial_scRNA_1,
                           reduction = 'pca',
                           group.by = 'orig.ident',
                           pt.size = 0.3) +
  ggtitle("Before Harmony: PCA by Sample") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

pdf('01_PCA_batch_before.pdf', width = 14, height = 6)
print(p_before + p_before_sample)
dev.off()

message("\n Harmony Info...")
Epithelial_scRNA_1 <- RunHarmony(Epithelial_scRNA_1,
                                 group.by.vars = "orig.ident",
                                 theta = 2,
                                 plot_convergence = FALSE)
message("Harmony Completed")

p_after <- DimPlot(Epithelial_scRNA_1,
                   reduction = 'harmony',
                   group.by = 'group',
                   pt.size = 0.5,
                   cols = c("BPH" = "#E64B35", "Healthy" = "#27AE60")) +
  ggtitle("After Harmony: Harmony Space by Group") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

p_after_sample <- DimPlot(Epithelial_scRNA_1,
                          reduction = 'harmony',
                          group.by = 'orig.ident',
                          pt.size = 0.3) +
  ggtitle("After Harmony: Harmony Space by Sample") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

pdf('02_PCA_batch_after.pdf', width = 14, height = 6)
print(p_after + p_after_sample)
dev.off()

pdf('03_PCA_batch_comparison.pdf', width = 14, height = 12)
print((p_before + p_after) / (p_before_sample + p_after_sample))
dev.off()

message("Info Generating: 01/02/03_PCA_batch_*.pdf")

message("\n Info (res 0.1-1.0)...")
plot_list <- list()
for (i in seq(0.1, 1, by = 0.1)) {
  Epithelial_temp <- Epithelial_scRNA_1 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = i) %>%
    RunUMAP(reduction = "harmony", dims = 1:10)

  p <- DimPlot(Epithelial_temp, reduction = "umap", pt.size = 0.3) +
    ggtitle(paste0("Resolution = ", i)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  plot_list[[as.character(i)]] <- p
}
rm(Epithelial_temp)

pdf('04_resolution_0.1_1.pdf', width = 28, height = 10)
wrap_plots(plot_list, ncol = 5)
dev.off()

message("\n Info res=0.6...")
Epithelial_scRNA_1 <- Epithelial_scRNA_1 %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = 0.6) %>%
  RunUMAP(reduction = "harmony", dims = 1:10)

pdf('05_Epithelial_umap_res0.6.pdf', width = 8, height = 6)
DimPlot(Epithelial_scRNA_1, label = T, pt.size = 0.5) +
  ggtitle("Epithelial Clusters (res=0.6, dims=1:10)") +
  theme_classic()
dev.off()

Epithelial_scRNA_1@meta.data$res_0.6 <- as.character(Idents(Epithelial_scRNA_1))

cols_to_keep <- c("orig.ident", "nCount_RNA", "nFeature_RNA",
                  "mt_percent", "group", "res_0.6")
existing_cols <- intersect(cols_to_keep, colnames(Epithelial_scRNA_1@meta.data))
Epithelial_scRNA_1@meta.data <- Epithelial_scRNA_1@meta.data[, existing_cols]

message("\n=== Computing Markers ===")
Epithelial_marker <- FindAllMarkers(Epithelial_scRNA_1,
                                    slot = 'data',
                                    only.pos = T,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)
write_csv(Epithelial_marker, 'Epithelial_marker.csv', col_names = T)

Epithelial_marker_top30 <- Epithelial_marker %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30, with_ties = FALSE) %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
write_csv(Epithelial_marker_top30, 'Epithelial_marker_top30.csv', col_names = T)

top30_wide <- Epithelial_marker %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30, with_ties = FALSE) %>%
  mutate(rank = row_number()) %>%
  select(cluster, gene, rank) %>%
  pivot_wider(names_from = rank, values_from = gene, names_prefix = "Top")
write_csv(top30_wide, 'Epithelial_marker_top30_wide.csv', col_names = T)

message("\n=== Generating Info genes Info ===")

epi_genes <- c('SCGB1A3', 'WFDC2', 'LCN2', 'MMP7', 'SCGB3A1',
               'KRT13', 'S100A16', 'S100A14', 'KRT19',
               'TP63', 'KRT14', 'KRT5',
               'KLK4', 'KLK3', 'KLK2', 'ACPP', 'AR')

pdf('06_Epithelial_genes_by_cluster_res0.6.pdf', width = 11, height = 8)
Idents(Epithelial_scRNA_1) <- Epithelial_scRNA_1$res_0.6

DotPlot(Epithelial_scRNA_1,
        features = epi_genes,
        group.by = "res_0.6") +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    name = "Avg Exp\n(scaled)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold", color = "black")
  ) +
  xlab('') + ylab('') +
  coord_flip()
dev.off()

message("\n=== Info ===")

Epithelial_scRNA_1 <- RenameIdents(Epithelial_scRNA_1,
                                   '0' = 'Basal',
                                   '1' = 'Basal',
                                   '2' = 'Basal',
                                   '3' = 'Hillock',
                                   '4' = 'Club',
                                   '5' = 'Hillock',
                                   '6' = 'Luminal_AR+',
                                   '7' = 'Basal',
                                   '8' = 'Luminal_AR-',
                                   '9' = 'Basal')

name <- as.character(unique(Idents(Epithelial_scRNA_1)))
if ("Unknown" %in% name) {
  name <- setdiff(name, "Unknown")
  Epithelial_scRNA_final <- Epithelial_scRNA_1[, Idents(Epithelial_scRNA_1) %in% name]
} else {
  Epithelial_scRNA_final <- Epithelial_scRNA_1
}

Epithelial_scRNA_final$cell_type <- factor(Idents(Epithelial_scRNA_final),
                                           levels = c("Club", "Hillock", "Basal",
                                                      "Luminal_AR+", "Luminal_AR-"))

Idents(Epithelial_scRNA_final) <- Epithelial_scRNA_final$cell_type

epithelial_colors <- c("Basal" = "#4DBBD5",
                       "Club" = "#91D1C2",
                       "Hillock" = "#E6B735",
                       "Luminal_AR+" = "#B2182B",
                       "Luminal_AR-" = "#F4A582")

pdf('07_Epithelial_anno_umap_final_res0.6.pdf', width = 10, height = 8)
DimPlot(Epithelial_scRNA_final,
        label = FALSE,
        label.size = 4,
        pt.size = 0.6,
        cols = epithelial_colors) +
  ggtitle("Epithelial Subtypes (res=0.6, Harmony Corrected, dims=1:10)") +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.8),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )
dev.off()

pdf('08_Epithelial_anno_dotplot_final_res0.6.pdf', width = 11, height = 9)
DotPlot(Epithelial_scRNA_final,
        features = unique(epi_genes),
        group.by = "cell_type") +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    name = "Avg Exp\n(scaled)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold", color = "black")
  ) +
  xlab('') + ylab('') +
  coord_flip()
dev.off()

message("\n cells Info:")
stats_table <- table(Epithelial_scRNA_final$cell_type, Epithelial_scRNA_final$group)
print(stats_table)

saveRDS(Epithelial_scRNA_final, 'Epithelial_bph_subcluster_res0.6.rds')
message("\n Completed! Generating file columns Info:")
message("01_PCA_batch_before.pdf - Info PCA")
message("02_PCA_batch_after.pdf - Info PCA")
message("03_PCA_batch_comparison.pdf - Info (4 Info)")
message("04_resolution_0.1_1.pdf - Info")
message("05_Epithelial_umap_res0.6.pdf - Info UMAP (res=0.6)")
message("06_Epithelial_genes_by_cluster_res0.6.pdf - Info genes Info")
message("07_Epithelial_anno_umap_final_res0.6.pdf - Info UMAP(None Info)")
message("08_Epithelial_anno_dotplot_final_res0.6.pdf - Info DotPlot")

output_dir <- "./scDRS_epithelial_res0.6"
if (!dir.exists(output_dir)) dir.create(output_dir)

library(Matrix)
DefaultAssay(Epithelial_scRNA_final) <- "RNA"
counts <- GetAssayData(Epithelial_scRNA_final, layer = "counts")

writeMM(counts, file.path(output_dir, "counts_raw.mtx"))
write.table(rownames(counts), file.path(output_dir, "genes.csv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(counts), file.path(output_dir, "cells.csv"), quote = FALSE, row.names = FALSE, col.names = FALSE)

meta <- data.frame(
  cell_id = colnames(counts),
  orig.ident = Epithelial_scRNA_final$orig.ident,
  group = Epithelial_scRNA_final$group,
  cell_type = as.character(Epithelial_scRNA_final$cell_type),
  nFeature_RNA = Epithelial_scRNA_final$nFeature_RNA,
  mt_percent = Epithelial_scRNA_final$mt_percent,
  stringsAsFactors = FALSE
)
write.csv(meta, file.path(output_dir, "metadata.csv"), quote = FALSE, row.names = FALSE)
message(sprintf("scDRS file Info Generating: %s/", output_dir))

sessionInfo()
