library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

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

output_dir <- "01_Epithelial_Subclustering"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Info cells Info (Info group Info)\n")
cat("========================================\n")

cat("\n[1/8] Info scRNA_harmony_annotated.rds...\n")
scRNA <- readRDS("scRNA_harmony_annotated.rds")

scRNA$group <- as.character(scRNA$group)
scRNA$group[scRNA$group == "Normal"] <- "Healthy"
scRNA$group <- ifelse(scRNA$group %in% c("Tumor_LowRisk", "Tumor_HighRisk"), "Tumor", scRNA$group)
scRNA$group <- factor(scRNA$group, levels = c("Healthy", "Tumor"))

cat("- Info cells Info:", ncol(scRNA), "\n")
cat("- Info group:\n")
print(table(scRNA$group, useNA = "ifany"))

cat("\n[2/8] Info cells...\n")
epi_all <- subset(scRNA, subset = cell_type == "Epithelial")

na_count <- sum(is.na(epi_all$group))
if (na_count > 0) {
  cat(sprintf("! Warning: group columns Info %.0f NA, Info \n", na_count))
  epi_all <- subset(epi_all, subset = !is.na(group))
}

cat("- Info cells:", ncol(epi_all), "\n")
cat("  - Healthy:", sum(epi_all$group == "Healthy", na.rm = TRUE), "\n")
cat("  - Tumor:", sum(epi_all$group == "Tumor", na.rm = TRUE), "\n")

cat("\n[3/8] rerunning Info and Info...\n")
epi_all <- epi_all %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

p_elbow <- ElbowPlot(epi_all, ndims = 30) +
  ggtitle("Epithelial Cells - Elbow Plot") +
  theme_clean()
ggsave(file.path(output_dir, "Figures", "01_ElbowPlot.pdf"),
       p_elbow, width = 6, height = 5, dpi = 300)

dims_use <- 1:10
cat("- Using PCs:", paste(dims_use, collapse = ","), "\n")

cat("\n[4/8] Info (k.param=20)...\n")
epi_all <- FindNeighbors(epi_all, dims = dims_use, k.param = 20, verbose = FALSE)

resolutions <- c(0.1, 0.15, 0.2, 0.3)
res_stats <- data.frame()

for (res in resolutions) {
  epi_all <- FindClusters(epi_all, resolution = res, verbose = FALSE)
  res_col <- paste0("RNA_snn_res.", res)
  n_clust <- length(unique(epi_all@meta.data[[res_col]]))
  clust_table <- table(epi_all@meta.data[[res_col]])

  res_stats <- rbind(res_stats, data.frame(
    Resolution = res,
    N_Clusters = n_clust,
    Min_Size = min(clust_table),
    Median_Size = median(as.numeric(clust_table)),
    Max_Size = max(clust_table)
  ))

  cat(sprintf("res=%.2f: %.0f Info | Info =%.0f | Info =%.0f | Info =%.0f cells \n",
              res, n_clust, min(clust_table), median(as.numeric(clust_table)), max(clust_table)))
}

print(res_stats)
write.csv(res_stats, file.path(output_dir, "02_Resolution_Stats.csv"), row.names = FALSE)

selected_res <- 0.3
selected_col <- paste0("RNA_snn_res.", selected_res)
epi_all$subcluster <- epi_all@meta.data[[selected_col]]

cat(sprintf("\n => Info resolution=%.2f, Info %.0f Info \n",
            selected_res, length(unique(epi_all$subcluster))))

cat("\n[5/8] Info rows UMAP...\n")
epi_all <- RunUMAP(epi_all, dims = dims_use,
                   n.neighbors = 12,
                   min.dist = 1.2,
                   spread = 2.0,
                   repulsion.strength = 2.0,
                   negative.sample.rate = 5,
                   verbose = FALSE)

cat("\n[6/8] Generating Info...\n")

group_colors_epi <- c("Healthy" = "#27AE60", "Tumor" = "#C0392B")
n_subclust <- length(unique(epi_all$subcluster))
cluster_colors <- rainbow(n_subclust, s = 0.8, v = 0.9)

p_cluster <- DimPlot(epi_all, group.by = "subcluster", label = TRUE,
                     repel = TRUE, pt.size = 0.5,
                     label.size = 5, label.box = TRUE, label.color = "black",
                     cols = cluster_colors) +
  ggtitle(sprintf("Epithelial Subclusters (res=%.2f)", selected_res)) +
  theme_clean() +
  theme(legend.position = "right", panel.grid = element_blank())
ggsave(file.path(output_dir, "Figures", "03_UMAP_Subcluster.pdf"),
       p_cluster, width = 10, height = 8, dpi = 300)

p_group <- DimPlot(epi_all, group.by = "group",
                   cols = group_colors_epi, pt.size = 0.5, label = FALSE) +
  ggtitle("Epithelial Cells by Group") +
  theme_clean() +
  theme(legend.position = "right")
ggsave(file.path(output_dir, "Figures", "04_UMAP_Group_TwoGroups.pdf"),
       p_group, width = 10, height = 8, dpi = 300)

p_split_group <- DimPlot(epi_all, group.by = "group",
                         split.by = "group",
                         cols = group_colors_epi,
                         pt.size = 0.4, ncol = 2) +
  ggtitle("Epithelial Cells - Split by Group") +
  theme_clean() +
  theme(strip.background = element_rect(fill = "black", color = NA),
        strip.text = element_text(size = 12, face = "bold", color = "white"),
        legend.position = "none")
ggsave(file.path(output_dir, "Figures", "05_UMAP_Split_by_Group.pdf"),
       p_split_group, width = 12, height = 6, dpi = 300)

p_split_sub <- DimPlot(epi_all, group.by = "group",
                       split.by = "subcluster",
                       cols = group_colors_epi,
                       pt.size = 0.5, ncol = 4) +
  ggtitle("Subcluster Composition (Two Groups)") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())
ggsave(file.path(output_dir, "Figures", "06_UMAP_Split_by_Subcluster.pdf"),
       p_split_sub, width = 16,
       height = 4 * ceiling(n_subclust/4), dpi = 300)

if(sum(epi_all$group == "Healthy", na.rm = TRUE) > 0) {
  p_healthy <- DimPlot(subset(epi_all, group == "Healthy"),
                       group.by = "subcluster",
                       label = TRUE, label.box = TRUE, pt.size = 1.5,
                       cols = cluster_colors) +
    ggtitle("Healthy Epithelial Cells Only") +
    theme_clean()
  ggsave(file.path(output_dir, "Figures", "07_UMAP_Healthy_Only.pdf"),
         p_healthy, width = 8, height = 6, dpi = 300)
}

p_combined <- DimPlot(epi_all, group.by = "subcluster", label = TRUE,
                      label.box = TRUE, pt.size = 0.3, cols = cluster_colors) +
  DimPlot(epi_all, group.by = "group",
          cols = group_colors_epi,
          pt.size = 0.3, alpha = 0.5) +
  plot_layout(ncol = 2) +
  plot_annotation(title = "Combined View: Subclusters & Two Groups")
ggsave(file.path(output_dir, "Figures", "08_UMAP_Combined.pdf"),
       p_combined, width = 14, height = 6, dpi = 300)

cat("\n[7/8] Info group Info (Info group)...\n")

comp_table2 <- table(epi_all$subcluster, epi_all$group)
comp_df2 <- as.data.frame.matrix(comp_table2)
comp_df2$Total <- rowSums(comp_df2)
comp_df2$Healthy_Pct <- round(comp_df2$Healthy / comp_df2$Total * 100, 1)
comp_df2$Tumor_Pct <- round(comp_df2$Tumor / comp_df2$Total * 100, 1)

comp_df2$Type <- apply(comp_df2[, c("Healthy_Pct", "Tumor_Pct")], 1, function(x) {
  if (x["Healthy_Pct"] >= 70) return("Healthy-rich")
  if (x["Tumor_Pct"] >= 70) return("Tumor-rich")
  return("Mixed")
})

cat("Info group group Info:\n")
print(comp_df2)
write.csv(comp_df2, file.path(output_dir, "09_Subcluster_Composition_TwoGroups.csv"))

comp_data2 <- as.data.frame(comp_table2)
colnames(comp_data2) <- c("Subcluster", "Group", "Count")

p_comp2 <- ggplot(comp_data2, aes(x = Subcluster, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = group_colors_epi) +
  geom_text(aes(label = ifelse(Count > 0, sprintf("%.0f", Count), "")),
            position = position_fill(vjust = 0.5), size = 3) +
  labs(title = "Subcluster Composition (Two Groups)",
       subtitle = sprintf("Total: %.0f cells", sum(comp_df2$Total)),
       x = "Subcluster", y = "Proportion") +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
ggsave(file.path(output_dir, "Figures", "10_Composition_Bar_TwoGroups.pdf"),
       p_comp2, width = 10, height = 6, dpi = 300)

cat("\n[8/8] Saved results...\n")
saveRDS(epi_all, file.path(output_dir, "13_epithelial_subclustered.rds"))
cat(sprintf("- results Info Saved Info: %s\n", output_dir))

cat("\n========================================\n")
cat("Info:\n")
cat(sprintf("- Resolution: %.2f (%.0f Info)\n", selected_res, n_subclust))
cat("Info cells Info:\n")
print(table(epi_all$subcluster))
cat("Info group group Info:\n")
print(comp_df2[, c("Healthy", "Tumor", "Total", "Type")])
cat("\n========================================\n")
cat("Info cells Info (Info group Info) Completed!\n")
cat("========================================\n")

sessionInfo()
