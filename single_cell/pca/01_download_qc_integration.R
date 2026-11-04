library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(R.utils)
library(data.table)
library(patchwork)
library(harmony)
library(SingleCellExperiment)
library(scDblFinder)
library(scales)
library(clustree)

rm(list=ls())
gc()

dir='GSE181294_RAW/'
samples <- list.files(dir)
print("Info Info:")
print(samples)

scRNAlist = lapply(samples, function(pro){
print(pro)
ct = fread(file.path(dir, pro), data.table = F)
rownames(ct) = ct[,1]
ct = ct[,-1]

sce = CreateSeuratObject(counts = ct,
                     project = gsub('.csv.gz$', '', pro),
                     min.cells = 5
)
return(sce)
})

for (i in seq_along(scRNAlist)) {
DefaultAssay(scRNAlist[[i]]) <- "RNA"
}

scRNAlist <- lapply(scRNAlist, function(x) {

mt_genes <- grep("^MT-|^MT", rownames(x), value = TRUE, ignore.case = TRUE)
if (length(mt_genes) > 0) {
x[["mt_percent"]] <- PercentageFeatureSet(x, features = mt_genes)
} else {
x[["mt_percent"]] <- 0
}

hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
hb_genes_found <- rownames(x)[toupper(rownames(x)) %in% toupper(hb_genes)]

if (length(hb_genes_found) > 0) {
x[["HB_percent"]] <- PercentageFeatureSet(x, features = hb_genes_found)
} else {
x[["HB_percent"]] <- runif(ncol(x), min = -0.001, max = 0.001)
}
return(x)
})

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

if(length(scRNAlist) == 1) {
scRNA <- scRNAlist[[1]]
} else {
scRNA <- Reduce(function(x, y) merge(x, y), scRNAlist)
}
message(sprintf("Merging Completed, Info cells Info: %d", ncol(scRNA)))

scRNA <- JoinLayers(scRNA)
message(sprintf("JoinLayers Completed"))

message("Info cells Info (scDblFinder, dbr=0.08)...")
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

n_samples <- length(unique(scRNA$orig.ident))
if(n_samples > 0) {
if(n_samples == 1) {
color_palette <- "#1f77b4"
} else {
color_palette <- scales::hue_pal()(n_samples)
}

meta_df <- as.data.frame(scRNA@meta.data)

p1 <- meta_df %>%
dplyr::count(orig.ident) %>%
ggplot(aes(x = reorder(orig.ident, -n), y = n, fill = orig.ident)) +
geom_col(show.legend = FALSE, width = 0.7) +
scale_fill_manual(values = color_palette) +
labs(title = "Cells per Sample (Doublets Removed, 8%)", x = "", y = "Number of Cells") +
theme_minimal() +
theme(
axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 1),
plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
)

p2 <- meta_df %>%
group_by(orig.ident) %>%
summarise(across(c(nFeature_RNA, nCount_RNA, mt_percent, HB_percent), median), .groups = "drop") %>%
tidyr::pivot_longer(-orig.ident, names_to = "Metric", values_to = "Value") %>%
ggplot(aes(x = orig.ident, y = Value, color = orig.ident)) +
geom_point(size = 3) +
facet_wrap(~Metric, scales = "free_y", ncol = 2) +
scale_color_manual(values = color_palette) +
labs(title = "Median Values per Sample", x = "", y = "Median Value") +
theme_minimal() +
theme(
axis.text.x = element_text(angle = 60, hjust = 1),
legend.position = "none",
plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
)

p3 <- meta_df %>%
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
plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
)

combined_plot <- p1 / p2 / p3
print(combined_plot)
}
dev.off()

message("Info and Info...")

scRNA <- NormalizeData(scRNA) %>%
FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
ScaleData() %>%
RunPCA(npcs = 30, verbose = T)

a <- DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")

top15 <- head(VariableFeatures(scRNA), 15)
plot_base <- VariableFeaturePlot(scRNA)
plot_labeled <- LabelPoints(plot = plot_base, points = top15, repel = TRUE, size = 3)
print(plot_labeled)
ggsave("top15_highly_variable_genes.pdf", plot = plot_labeled, width = 12, height = 8)

g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))

scRNA <- CellCycleScoring(object = scRNA,
                    g2m.features = g2m_genes,
                    s.features = s_genes)

scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident", plot_convergence = T)

harmony_plot <- DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident")
pca_harmony_integrated <- a + harmony_plot + plot_layout(ncol = 2) +
plot_annotation(title = "Batch Effect Removal Comparison")
print(pca_harmony_integrated)
ggsave("pca_harmony_comparison.pdf", pca_harmony_integrated, width = 14, height = 7)

elbow_plot <- ElbowPlot(scRNA_harmony, ndims = 30, reduction = "harmony") +
geom_vline(xintercept = 10, linetype = "dashed", color = "red")
print(elbow_plot)
ggsave("elbow_plot.pdf", elbow_plot, width = 8, height = 6)

scRNA_harmony <- FindNeighbors(scRNA_harmony, dims = 1:10, reduction = "harmony")

resolutions <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0)
for (res in resolutions) {
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = res)
}

clustree_plot <- clustree(scRNA_harmony, prefix = "RNA_snn_res.") +
ggtitle("Cluster Relationships")
print(clustree_plot)
ggsave("clustree_plot.pdf", clustree_plot, width = 10, height = 10)

scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.8, algorithm = 1)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:10)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:10)

umap_cluster <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
umap_batch <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")

combined_plot <- umap_cluster + umap_batch
print(combined_plot)
ggsave("umap_tsne_integrated.pdf", combined_plot, width = 14, height = 6)

cluster_counts <- table(scRNA_harmony@meta.data$seurat_clusters)
print("Info:")
print(cluster_counts)

save(scRNA_harmony, scRNAlist, file = "scdata2.Rdata")

message(sprintf("Completed! Total %d cells, %d Info",
          ncol(scRNA_harmony), length(unique(Idents(scRNA_harmony)))))

sessionInfo()
