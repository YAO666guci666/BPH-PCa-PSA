library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(R.utils)
library(data.table)
library(patchwork)
library(hdf5r)
library(harmony)
library(SingleCellExperiment)
library(scDblFinder)
library(scales)
library(clustree)

rm(list=ls())

gc()

getwd()

.libPaths()

dir='GSM5252137_RAW/'
samples <- c("BPH283PrGF","BPH327PrGF","BPH340PrGF","BPH389PrGF",
             "BPH511PrG","BPH556PrGA1","BPH556PrGA2",
             "D17PrTzF","D27PrTzF","D35PrTzF","D17PrPzF","D27PrPzF","D35PrPzF")

scRNAlist = lapply(samples, function(pro){
  print(pro)
  h5_file <- file.path(dir, paste0(pro, ".h5"))
  ct <- Read10X_h5(h5_file)
  colnames(ct) <- paste0(pro, "_", colnames(ct))

  sce = CreateSeuratObject(counts = ct,
                           project = pro,
                           min.cells = 3
  )

  sce$group <- ifelse(grepl("^BPH", pro), "BPH", "Healthy")

  return(sce)
})
scRNAlist
samples

for (i in seq_along(scRNAlist)) {
  DefaultAssay(scRNAlist[[i]]) <- "RNA"
}

scRNAlist <- lapply(scRNAlist, function(x) {

  mt_genes <- grep("^MT-", rownames(x), value = TRUE, ignore.case = TRUE)
  if (length(mt_genes) > 0) {
    x[["mt_percent"]] <- PercentageFeatureSet(x, features = mt_genes)
  } else {
    warning(paste("No mitochondrial genes found in sample:", x@project.name))
    x[["mt_percent"]] <- 0
  }

  hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

  hb_genes_found <- rownames(x)[toupper(rownames(x)) %in% toupper(hb_genes)]

  if (length(hb_genes_found) > 0) {
    x[["HB_percent"]] <- PercentageFeatureSet(x, features = hb_genes_found)
  } else {
    warning(paste("No hemoglobin genes found in sample:", x@project.name))

    x[["HB_percent"]] <- runif(ncol(x), min = -0.001, max = 0.001)
  }
  return(x)
})

pdf("qc_violins_before_filtering.pdf", width = 10, height = 6)

uniform_theme <- theme(
  axis.text = element_text(size = 7),
  axis.title = element_text(size = 8),
  plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
  legend.position = "none",
  plot.margin = unit(c(2, 2, 2, 2), "mm")
)

feature_labels <- c(
  "nFeature_RNA" = "Genes per cell",
  "nCount_RNA" = "UMIs per cell",
  "mt_percent" = "MT gene %",
  "HB_percent" = "HB gene %"
)

for (i in seq_along(scRNAlist)) {

  sample_name <- scRNAlist[[i]]@project.name
  cat("Plotting sample:", sample_name, "\n")

  features <- c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent")

  p_list <- VlnPlot(
    scRNAlist[[i]],
    features = features,
    pt.size = 0.1,
    combine = FALSE,
    fill.by = "feature"
  )

  for (j in seq_along(p_list)) {
    feature_name <- features[j]

    if (feature_name == "HB_percent") {
      hb_values <- scRNAlist[[i]]$HB_percent

      y_limits <- if (all(abs(hb_values) < 0.01)) {
        c(-0.1, 0.1)
      } else {
        range(hb_values)
      }

      p_list[[j]] <- p_list[[j]] +
        uniform_theme +
        ggtitle(feature_labels[feature_name]) +
        ylab("Value") +
        theme(axis.title.x = element_blank()) +
        scale_fill_manual(values = "#1f77b4") +
        guides(fill = "none") +
        ylim(y_limits)
    } else {
      p_list[[j]] <- p_list[[j]] +
        uniform_theme +
        ggtitle(feature_labels[feature_name]) +
        ylab("Value") +
        theme(axis.title.x = element_blank()) +
        scale_fill_manual(values = "#1f77b4") +
        guides(fill = "none")
    }
  }

  combined <- wrap_plots(p_list, ncol = 4) +
    plot_annotation(
      title = sample_name,
      theme = theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold"))
    )

  print(combined)
}
dev.off()

scRNAlist <- lapply(scRNAlist, function(x) {
  upper_quantile <- quantile(x$nCount_RNA, 0.97)
  subset(x,
         nFeature_RNA > 200 & nFeature_RNA < 5000 &
           mt_percent < 15 &
           HB_percent < 3 &
           nCount_RNA < upper_quantile &
           nCount_RNA > 1000)
})

message("Merging Info...")
scRNA <- merge(x = scRNAlist[[1]], y = scRNAlist[-1])

scRNA <- JoinLayers(scRNA)
message(sprintf("Merging Completed, Info cells Info: %d", ncol(scRNA)))

message("Info cells Info (scDblFinder)...")
DefaultAssay(scRNA) <- "RNA"
counts_mat <- GetAssayData(scRNA, layer = "counts")

sce <- SingleCellExperiment(
  assays = list(counts = counts_mat),
  colData = scRNA@meta.data
)

sce <- scDblFinder(sce, samples = "orig.ident", dbr = 0.08, nfeatures = 1000)
scRNA$scDblFinder.class <- sce$scDblFinder.class
scRNA$scDblFinder.score <- sce$scDblFinder.score

rm(sce, counts_mat); gc()

before <- ncol(scRNA)
scRNA <- subset(scRNA, scDblFinder.class == "singlet")
after <- ncol(scRNA)
message(sprintf("Info cells Info: %d > %d (Info %.1f%%)", before, after, 100*(before-after)/before))

pdf("qc_summary_after_filtering.pdf", width = 16, height = 14)

color_palette <- scales::hue_pal()(length(unique(scRNA$orig.ident)))
num_samples <- length(unique(scRNA$orig.ident))

p1 <- scRNA@meta.data %>%
  as.data.frame() %>%
  dplyr::count(orig.ident) %>%
  ggplot(aes(x = reorder(orig.ident, -n), y = n, fill = orig.ident)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  scale_fill_manual(values = color_palette) +
  labs(title = "Cells per Sample After QC & Doublet Removal", x = "", y = "Number of Cells") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 1),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.margin = margin(t = 10, r = 20, b = 20, l = 10)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

p2 <- scRNA@meta.data %>%
  as.data.frame() %>%
  group_by(orig.ident) %>%
  summarise(across(c(nFeature_RNA, nCount_RNA, mt_percent, HB_percent), median), .groups = "drop") %>%
  tidyr::pivot_longer(-orig.ident, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = orig.ident, y = Value, color = orig.ident)) +
  geom_point(size = 3) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_color_manual(values = color_palette) +
  scale_y_continuous(limits = ~ c(min(.x, na.rm = TRUE),
                                  max(.x, na.rm = TRUE) * 1.25)) +
  labs(title = "Median Values per Sample", x = "", y = "Median Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10, vjust = 1),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1.5, "lines")
  )

p3 <- scRNA@meta.data %>%
  as.data.frame() %>%
  tidyr::pivot_longer(c(nFeature_RNA, nCount_RNA, mt_percent, HB_percent),
                      names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Value, color = orig.ident)) +
  geom_density(alpha = 0.7, linewidth = 0.8) +
  facet_wrap(~Metric, scales = "free", ncol = 2) +
  scale_color_manual(values = color_palette) +
  labs(title = "Distribution of QC Metrics", x = "Value", y = "Density") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = ifelse(num_samples > 10, 7, 9)),
    legend.key.size = unit(ifelse(num_samples > 10, 0.6, 0.8), "lines"),
    legend.box.spacing = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
  ) +
  guides(color = guide_legend(
    nrow = ifelse(num_samples > 10, 3, ifelse(num_samples > 5, 2, 1)),
    byrow = TRUE,
    override.aes = list(size = ifelse(num_samples > 10, 2, 3), alpha = 1)
  ))

combined_plot <- p1 / p2 / p3 +
  plot_layout(heights = c(0.8, 1.5, ifelse(num_samples > 10, 3, 2.5))) &
  theme(plot.margin = margin(20, 20, ifelse(num_samples > 10, 40, 30), 20))

print(combined_plot)
dev.off()

scRNAlist <- NormalizeData(scRNA) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = T)

a <- DimPlot(scRNAlist, reduction = "pca", group.by = "orig.ident")
a

top15 <- head(VariableFeatures(scRNAlist), 15)

plot_base <- VariableFeaturePlot(scRNAlist)

plot_labeled <- LabelPoints(
  plot = plot_base,
  points = top15,
  repel = TRUE,
  size = 3,
  xnudge = 0.1,
  ynudge = 0.1
)

print(plot_labeled)

ggsave(
  filename = "top15_highly_variable_genes.pdf",
  plot = plot_labeled,
  width = 12,
  height = 8
)

scRNAlist <- JoinLayers(scRNAlist)

g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNAlist))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNAlist))

scRNAlist <- CellCycleScoring(object = scRNAlist,
                              g2m.features = g2m_genes,
                              s.features = s_genes,
                              set.ident = TRUE)

p4 <- VlnPlot(scRNAlist, features = c("S.Score", "G2M.Score"),
              group.by = "orig.ident", ncol = 2, pt.size = 0.1)
p4

scRNAlist@meta.data %>%
  as.data.frame() %>%
  ggplot(aes(S.Score, G2M.Score)) +
  geom_point(aes(color = Phase)) +
  theme_minimal()

scRNA_harmony <- RunHarmony(
  scRNAlist,
  group.by.vars = "orig.ident",
  plot_convergence = TRUE
)

print(scRNA_harmony@reductions[["harmony"]]@cell.embeddings[1:5, 1:5])

harmony_plot <- DimPlot(
  scRNA_harmony,
  reduction = "harmony",
  group.by = "orig.ident",
  pt.size = 0.5
) +
  ggtitle("Harmony Integration - Batch Effect Removal")

pca_harmony_integrated <- a + harmony_plot +
  plot_layout(ncol = 2) +
  plot_annotation(title = "Batch Effect Removal Comparison")
print(pca_harmony_integrated)
ggsave("pca_harmony_comparison.pdf", pca_harmony_integrated, width = 14, height = 7)

elbow_plot <- ElbowPlot(scRNA_harmony, ndims = 30, reduction = "harmony") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  annotate("text", x = 12, y = 5, label = "Selected dims: 1-10", color = "red")
print(elbow_plot)
ggsave("elbow_plot.pdf", elbow_plot, width = 8, height = 6)

scRNA_harmony <- FindNeighbors(scRNA_harmony, dims = 1:10, reduction = "harmony")

resolutions <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0)
for (res in resolutions) {
  scRNA_harmony <- FindClusters(scRNA_harmony, resolution = res)
}

clustree_plot <- clustree(scRNA_harmony, prefix = "RNA_snn_res.") +
  ggtitle("Cluster Relationships") +
  theme(plot.title = element_text(hjust = 0.5, size = 25))

print(clustree_plot)
ggsave("clustree_plot.pdf", clustree_plot, width = 10, height = 10)

my_res <- 0.2
scRNA_harmony <- FindClusters(
  scRNA_harmony,
  resolution = my_res,
  algorithm = 1
)

scRNA_harmony <- RunTSNE(
  scRNA_harmony,
  reduction = "harmony",
  dims = 1:10,
  perplexity = 30
)
scRNA_harmony <- RunUMAP(
  scRNA_harmony,
  reduction = "harmony",
  dims = 1:10,
  n.neighbors = 30
)

create_dim_plot <- function(reduction, group.by = NULL, label = FALSE, cols = NULL) {
  if (is.null(group.by)) group.by <- "seurat_clusters"
  p <- DimPlot(scRNA_harmony,
               reduction = reduction,
               group.by = group.by,
               label = label,
               cols = cols) +
    theme_minimal() +
    labs(title = paste0(toupper(reduction), " - ",
                        ifelse(label, "Clusters", "Batch")))
  if (label) p <- p + NoLegend()
  p
}

tsne_batch   <- create_dim_plot("tsne", group.by = "orig.ident")
tsne_cluster <- create_dim_plot("tsne", label = TRUE)
umap_batch   <- create_dim_plot("umap", group.by = "orig.ident")
umap_cluster <- create_dim_plot("umap", label = TRUE)

umap_group <- create_dim_plot("umap", group.by = "group",
                              cols = c("BPH"="#E64B35", "Healthy"="#00A087")) +
  ggtitle("UMAP - Group (BPH vs Healthy)")

umap_tsne_integrated <- (tsne_batch + tsne_cluster) /
  (umap_batch + umap_cluster) +
  plot_annotation(title = "Integration Results",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

print(umap_tsne_integrated)
ggsave("umap_tsne_integrated.pdf", umap_tsne_integrated, width = 14, height = 12)

ggsave("umap_group.pdf", umap_group, width = 8, height = 7)

cluster_counts <- table(scRNA_harmony@meta.data$seurat_clusters)
print("Cluster cell counts (resolution 0.2):")
print(cluster_counts)

save(scRNA_harmony, scRNAlist, file = "scdata2.Rdata")

cluster_bar <- ggplot(data.frame(cluster_counts),
                      aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq), vjust = -0.3) +
  labs(x = "Cluster", y = "Cell Count",
       title = "Cell Distribution per Cluster (res 0.2)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")
print(cluster_bar)
ggsave("cluster_distribution.pdf", cluster_bar, width = 10, height = 6)

message("data Info, Info and Info cells Info Completed!")
message(sprintf("Info cells Info: %d, genes Info: %d", ncol(scRNA_harmony), nrow(scRNA_harmony)))

sessionInfo()
