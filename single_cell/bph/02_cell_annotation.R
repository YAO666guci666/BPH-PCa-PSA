options(stringsAsFactors = FALSE)

rm(list=ls())

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(scales)
load("./scdata2.Rdata")

plot1 = DimPlot(scRNA_harmony, reduction = "umap",label = T)
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident')
plot1 = DimPlot(scRNA_harmony, reduction = "tsne",label = T)
plot2 = DimPlot(scRNA_harmony, reduction = "tsne", group.by='orig.ident')

scRNA_harmony.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE,
                                        min.pct = 0.25, logfc.threshold = 0.25)

top30 <- scRNA_harmony.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)

write.csv(top30, "cluster_30markers.csv", row.names = FALSE)

theme_clean <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.line = element_line(color = "black", size = 0.8),
      axis.text = element_text(color = "black", size = 14),
      axis.title = element_text(color = "black", size = 16, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black", size = 0.8)
    )
}

pdf("01_UMAP_Clusters_Before_Manual_Annotation.pdf", width = 12, height = 10)
cluster_colors <- rainbow(length(unique(scRNA_harmony@meta.data$seurat_clusters)),
                          s = 0.8, v = 0.9, alpha = 0.8)
names(cluster_colors) <- sort(unique(scRNA_harmony@meta.data$seurat_clusters))

DimPlot(scRNA_harmony, reduction = "umap",
        pt.size = 0.8,
        label = TRUE,
        label.size = 6,
        label.color = "black",
        label.box = TRUE,
        cols = cluster_colors) +
  ggtitle("UMAP - Cell Clusters (Before Manual Annotation)") +
  theme_clean() +
  theme(legend.position = "right") +
  guides(color = guide_legend(title = "Clusters",
                              override.aes = list(size = 4),
                              ncol = 1))
dev.off()

message("Info UMAP Info Generating: 01_UMAP_Clusters_Before_Manual_Annotation.pdf")

top5_clean <- scRNA_harmony.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  distinct(gene, .keep_all = TRUE)

DefaultAssay(scRNA_harmony) <- "RNA"
exp_data <- FetchData(scRNA_harmony, vars = unique(top5_clean$gene))
exp_data$cluster <- Idents(scRNA_harmony)
if ("seurat_clusters" %in% colnames(scRNA_harmony@meta.data)) {
  exp_data$cluster <- scRNA_harmony$seurat_clusters
}

avg_exp <- exp_data %>%
  group_by(cluster) %>%
  summarise(across(all_of(unique(top5_clean$gene)), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = -cluster, names_to = "gene", values_to = "avg_exp")

pct_exp <- exp_data %>%
  group_by(cluster) %>%
  summarise(across(all_of(unique(top5_clean$gene)), ~ mean(.x > 0, na.rm = TRUE))) %>%
  pivot_longer(cols = -cluster, names_to = "gene", values_to = "pct_exp")

plot_data <- left_join(avg_exp, pct_exp, by = c("cluster", "gene"))

gene_order <- top5_clean %>% arrange(cluster) %>% pull(gene)
plot_data$gene <- factor(plot_data$gene, levels = rev(gene_order))
plot_data$cluster <- factor(plot_data$cluster)

p <- ggplot(plot_data, aes(x = cluster, y = gene, size = pct_exp, color = avg_exp)) +
  geom_point() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  scale_size(range = c(1, 8)) +
  theme_classic(base_size = 12) +
  labs(title = "Top 5 Marker Genes per Cluster (DotPlot)",
       x = "Cluster",
       y = "Gene",
       color = "Average\nExpression",
       size = "Percent\nExpressed") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

pdf("02_Cluster_Markers_DotPlot.pdf", width = 14, height = 8)
print(p)
dev.off()

get_expected_markers <- function() {
  validation_genes <- list(
    "T Cells" = c("CD3D", "CD3E", "CD3G"),
    "Smooth Muscle" = c("ACTA2", "MYH11", "NDUFA4L2", "RGS5"),
    "Endothelial" = c("AQP1", "EMCN", "PLVAP", "VWF"),
    "Myeloid Cells" = c("CD14", "CD68", "LYZ"),
    "Fibroblasts" = c("APOD", "COL1A2", "DCN", "LUM"),
    "Mast Cells" = c("CPA3", "IL1RL1", "KIT", "TPSAB1"),
    "B Cells" = c("BANK1", "CD79A", "CD79B", "MS4A1"),
    "Epithelial" = c("EPCAM", "KRT8", "KRT18", "CDH1")
  )

  marker_list <- list()
  for (cell_type in names(validation_genes)) {
    genes <- validation_genes[[cell_type]]
    for (gene in genes) {
      marker_list[[length(marker_list) + 1]] <- data.frame(
        Cell_Type = cell_type,
        genes = gene
      )
    }
  }

  markers <- bind_rows(marker_list) %>%
    dplyr::filter(genes %in% rownames(scRNA_harmony)) %>%
    group_by(Cell_Type) %>%
    slice_head(n = 7) %>%
    distinct() %>%
    group_by(genes) %>%
    slice_head(n = 1) %>%
    ungroup()

  return(markers)
}

generate_expected_dot_plot <- function(marker_data, seurat_obj) {

  gene_list <- unique(marker_data$genes)

  dot_data <- DotPlot(seurat_obj,
                      features = gene_list,
                      group.by = "seurat_clusters",
                      dot.scale = 6)$data

  dot_data <- dot_data %>%
    left_join(marker_data, by = c("features.plot" = "genes"))

  gene_order <- dot_data %>%
    arrange(Cell_Type, features.plot) %>%
    distinct(features.plot, .keep_all = TRUE) %>%
    mutate(features.plot = factor(features.plot, levels = unique(features.plot)))

  dot_data$features.plot <- factor(dot_data$features.plot, levels = levels(gene_order$features.plot))

  cluster_order <- sort(unique(as.character(dot_data$id)))
  dot_data$id <- factor(dot_data$id, levels = cluster_order)

  cell_types <- unique(marker_data$Cell_Type)
  n_colors <- length(cell_types)

  saturated_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#999999", "#66C2A5",
    "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"
  )

  if (n_colors > length(saturated_colors)) {
    saturated_colors <- colorRampPalette(saturated_colors)(n_colors)
  } else {
    saturated_colors <- saturated_colors[1:n_colors]
  }

  names(saturated_colors) <- cell_types

  dot_colors <- colorRampPalette(c("#4393C3", "#F7F7F7", "#D6604D"))(100)

  expected_plot <- ggplot(dot_data, aes(x = id, y = features.plot)) +
    geom_point(aes(color = avg.exp.scaled, size = pct.exp),
               shape = 19, alpha = 0.8) +

    scale_color_gradientn(colors = dot_colors,
                          name = "Scaled\nExpression",
                          breaks = c(-2, 0, 2),
                          labels = c("Low", "Medium", "High")) +
    scale_size_continuous(range = c(1, 8),
                          breaks = c(0, 25, 50, 75, 100),
                          labels = c("0%", "25%", "50%", "75%", "100%"),
                          name = "Expression\nPercentage") +

    labs(title = "Expected Marker Genes Expression by Cluster",
         subtitle = "For cell-type annotation",
         x = "Cluster",
         y = "Gene") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(
        size = 10,
        color = saturated_colors[as.character(dot_data$Cell_Type[match(levels(dot_data$features.plot), dot_data$features.plot)])],
        face = "bold"
      ),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )

  pdf("03_Expected_Marker_DotPlot.pdf", width = 12, height = max(8, 0.4 * length(gene_list)))
  print(expected_plot)
  dev.off()

  return(expected_plot)
}

message("Info 1: Generating Info set genes Info...")
expected_markers <- get_expected_markers()

if(nrow(expected_markers) == 0) {
  stop("Not found Info Info set genes, check Info set genes Info")
}

expected_plot <- generate_expected_dot_plot(expected_markers, scRNA_harmony)
message("Info Generating: 03_Expected_Marker_DotPlot.pdf")

message("Info 2: Info rows cells Info...")

EXPECTED_CLUSTER_MAPPING <- c(
  "0" = "Epithelial",
  "1" = "Endothelial",
  "2" = "Epithelial",
  "3" = "T Cells",
  "4" = "B Cells",
  "5" = "Smooth Muscle",
  "6" = "Fibroblasts",
  "7" = "Epithelial",
  "8" = "Mast Cells",
  "9" = "Myeloid Cells",
  "10" = "Epithelial"
)

EXPECTED_CELLTYPE_ORDER <- c("Myeloid Cells", "Endothelial","Epithelial","Fibroblasts","Mast Cells","B Cells",
                             "Smooth Muscle" , "T Cells"
)

cluster_ids <- as.character(scRNA_harmony$seurat_clusters)
cell_type <- character(length(cluster_ids))

for (i in seq_along(cluster_ids)) {
  cluster_id <- cluster_ids[i]
  if (cluster_id %in% names(EXPECTED_CLUSTER_MAPPING)) {
    cell_type[i] <- EXPECTED_CLUSTER_MAPPING[[cluster_id]]
  } else {
    cell_type[i] <- paste0("Unknown_Cluster_", cluster_id)
    message(sprintf("Warning: Info %s Not found Info, Info set Info Unknown", cluster_id))
  }
}

scRNA_harmony <- AddMetaData(scRNA_harmony, metadata = cell_type, col.name = "cell_type")
message("cells Info Completed!")

cluster_celltype_mapping <- data.frame(
  Cluster = names(EXPECTED_CLUSTER_MAPPING),
  CellType = unlist(EXPECTED_CLUSTER_MAPPING),
  stringsAsFactors = FALSE
)

cluster_counts <- table(scRNA_harmony$seurat_clusters)
cluster_celltype_mapping$Cell_Count <- cluster_counts[cluster_celltype_mapping$Cluster]
cluster_celltype_mapping$Percentage <- round(
  cluster_celltype_mapping$Cell_Count / sum(cluster_counts) * 100, 2
)

write.csv(cluster_celltype_mapping, "Cluster_CellType_Mapping.csv", row.names = FALSE)

generate_annotated_dot_plot <- function(marker_data, seurat_obj,
                                        celltype_order = NULL) {

  gene_list <- unique(marker_data$genes)

  Idents(seurat_obj) <- seurat_obj$cell_type

  dot_data <- DotPlot(seurat_obj,
                      features = gene_list,
                      dot.scale = 6)$data

  dot_data <- dot_data %>%
    left_join(marker_data, by = c("features.plot" = "genes"))

  if(is.null(celltype_order)) {
    celltype_order <- EXPECTED_CELLTYPE_ORDER
  }

  existing_celltypes <- intersect(celltype_order, unique(dot_data$id))
  dot_data$id <- factor(dot_data$id, levels = existing_celltypes)

  gene_order <- dot_data %>%
    arrange(Cell_Type, features.plot) %>%
    distinct(features.plot, .keep_all = TRUE) %>%
    pull(features.plot)

  dot_data$features.plot <- factor(dot_data$features.plot, levels = gene_order)

  cell_types <- unique(marker_data$Cell_Type)
  n_colors <- length(cell_types)

  if (n_colors <= 8) {
    celltype_colors <- brewer.pal(8, "Set2")[1:n_colors]
  } else {
    celltype_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_colors)
  }

  names(celltype_colors) <- cell_types

  exp_colors <- colorRampPalette(c("#4575B4", "#F7F7F7", "#D73027"))(100)

  annotated_plot <- ggplot(dot_data, aes(x = id, y = features.plot)) +
    geom_point(aes(color = avg.exp.scaled, size = pct.exp),
               shape = 19, alpha = 0.8, stroke = 0.5) +

    scale_color_gradientn(
      colors = exp_colors,
      name = "Z-score",
      breaks = c(-2, 0, 2),
      labels = c("Low", "Medium", "High")
    ) +
    scale_size_continuous(
      range = c(2, 8),
      breaks = c(20, 40, 60, 80),
      labels = c("20%", "40%", "60%", "80%"),
      name = "Expression\nPercentage"
    ) +

    labs(
      title = "Cell Type Marker Expression Profile",
      subtitle = "Based on manual cell-type annotation",
      x = "Cell Type",
      y = "Marker Genes"
    ) +

    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5,
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey40",
                                   margin = margin(b = 15)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                 face = "bold", color = "black"),
      axis.text.y = element_text(
        size = 11,
        color = celltype_colors[as.character(dot_data$Cell_Type[match(levels(dot_data$features.plot),
                                                                      dot_data$features.plot)])],
        face = "bold"
      ),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(1, 1, 1, 2, "cm")
    ) +

    geom_vline(xintercept = seq(1.5, length(unique(dot_data$id)) - 0.5, 1),
               color = "grey80", size = 0.5, linetype = "dashed")

  pdf("04_Annotated_CellType_DotPlot.pdf", width = 14,
      height = max(8, 0.3 * length(gene_list)))
  print(annotated_plot)
  dev.off()

  return(annotated_plot)
}

message("Info 3: Generating Info cells Info...")
annotated_plot <- generate_annotated_dot_plot(expected_markers, scRNA_harmony,
                                              celltype_order = EXPECTED_CELLTYPE_ORDER)
message("Info cells Info Generating: 04_Annotated_CellType_DotPlot.pdf")

message("Info 4: Generating Info Info...")

n_manual_types <- length(unique(scRNA_harmony$cell_type))

if (n_manual_types <= 8) {
  manual_colors <- brewer.pal(8, "Set2")[1:n_manual_types]
} else if (n_manual_types <= 12) {
  manual_colors <- brewer.pal(12, "Set3")[1:n_manual_types]
} else {
  manual_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#98df8a", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
  )[1:n_manual_types]
}

names(manual_colors) <- sort(unique(scRNA_harmony$cell_type))

message("Info 5: Generating Info group UMAP Info...")

sample_ids <- scRNA_harmony$orig.ident

group_vec <- ifelse(grepl("^BPH", sample_ids, ignore.case = TRUE), "BPH",
                    "Healthy")

scRNA_harmony$group <- group_vec
group_col <- "group"

message("Info group Completed:")
message("  - BPH: ", sum(group_vec == "BPH"), "Info")
message("  - Healthy: ", sum(group_vec == "Healthy"), "Info")

group_colors <- c(
  "BPH" = "#E64B35",
  "Healthy" = "#27AE60"
)

pdf("07_UMAP_Group_Distribution.pdf", width = 14, height = 7)

p_group <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "group",
                   cols = group_colors, pt.size = 0.5, label = FALSE) +
  ggtitle("Group Distribution (BPH vs Healthy)") +
  theme_clean() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1))

p_split <- DimPlot(scRNA_harmony, reduction = "umap", split.by = "group",
                   group.by = "group", cols = group_colors,
                   pt.size = 0.4, ncol = 2) +
  ggtitle("Split View by Group") +
  theme_clean() +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 12, face = "bold", color = "white"),
        strip.background = element_rect(fill = "black", color = NA),
        legend.position = "none")

print(p_group + p_split)
dev.off()

pdf("08_UMAP_CellType_by_Group.pdf", width = 12, height = 5.5)

p1 <- DimPlot(scRNA_harmony, reduction = "umap", split.by = "group",
              group.by = "cell_type",
              cols = manual_colors,
              pt.size = 0.3,
              ncol = 2,
              raster=FALSE,
              label = TRUE,
              label.size = 2.5,
              repel = TRUE) +
  ggtitle("Cell Type Distribution: BPH vs Healthy") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(fill = "white", color = "black", size = 0.8),
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    panel.spacing = unit(1, "lines")
  ) +
  NoLegend()

print(p1)
dev.off()

message("Info group UMAP Info Generating: 07_UMAP_Group_Distribution.pdf, 08_UMAP_CellType_by_Group.pdf")
message("Info group Info:")
print(table(scRNA_harmony$group, useNA = "ifany"))

saveRDS(scRNA_harmony, "scRNA_harmony_annotated.rds")
message("Info Seurat Info Saved: scRNA_harmony_annotated.rds")

message("\n Info cells Info:")
print(table(scRNA_harmony$cell_type))
message("cells Info and Info group Info Completed!")

library(Matrix)
DefaultAssay(scRNA_harmony) <- "RNA"
counts <- GetAssayData(scRNA_harmony, layer = "counts")

writeMM(counts, "counts_raw.mtx")
write.table(rownames(counts), "genes.csv", quote=F, row.names=F, col.names=F)
write.table(colnames(counts), "cells.csv", quote=F, row.names=F, col.names=F)

meta_cols <- c("orig.ident", "seurat_clusters", "cell_type", "group",
               "nCount_RNA", "nFeature_RNA", "mt_percent")
meta_cols <- intersect(meta_cols, colnames(scRNA_harmony@meta.data))

meta <- data.frame(
  cell_id = colnames(counts),
  scRNA_harmony@meta.data[, meta_cols, drop = FALSE]
)
write.csv(meta, "metadata.csv", quote=F, row.names=F)

message("scDRS file Info Generating: counts_raw.mtx, genes.csv, cells.csv, metadata.csv")

sessionInfo()
