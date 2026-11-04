rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(dplyr)
library(ggplot2)

OUT_DIR <- "03_Epithelial_Annotation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Info data...\n")
epi <- readRDS("02_CopyKAT_Combined_Tumor/epi_all_copykat_annotated.rds")

if (any(duplicated(colnames(epi)))) {
  colnames(epi) <- make.names(colnames(epi), unique = TRUE)
}
rownames(epi@meta.data) <- colnames(epi)

for (a in names(epi@assays)) {
  for (s in c("counts", "data", "scale.data")) {
    if (s %in% slotNames(epi@assays[[a]])) {
      m <- slot(epi@assays[[a]], s)
      if (!is.null(m) && ncol(m) > 0) colnames(m) <- colnames(epi)
      slot(epi@assays[[a]], s) <- m
    }
  }
}

epi <- epi[, !is.na(epi$malignant_status)]
cat("Info cells:", ncol(epi), "\n")

if (!"data" %in% names(epi@assays$RNA)) epi <- NormalizeData(epi, verbose = FALSE)
epi$subcluster <- as.character(epi$subcluster)

epi_markers <- list(
  Basal   = c("KRT5", "KRT14", "TP63", "NGFR"),
  Luminal = c("KLK2", "KLK3", "KLK4", "ACPP"),
  Club    = c("SCGB3A1", "SCGB1A1", "MMP7", "LCN2")
)
all_markers <- unlist(epi_markers, use.names = FALSE)
all_markers <- all_markers[all_markers %in% rownames(epi)]

cat("Using marker:\n")
print(epi_markers)

cat("\n Plotting Info...\n")

normal_cells <- colnames(epi)[
  epi$malignant_status %in% c("Normal-like", "Reference") &
    !(epi$subcluster %in% c("5", "6"))
]

epi_normal <- subset(epi, cells = normal_cells)
epi_normal$subcluster <- factor(paste0("C", epi_normal$subcluster))
Idents(epi_normal) <- "subcluster"

p1 <- DotPlot(epi_normal, features = all_markers, cols = c("white", "#E64B35"), dot.scale = 7) +
  coord_flip() +
  labs(
    title = "Marker Expression by Subcluster (Before Annotation)",
    subtitle = "Normal subclusters only (C5/C6 excluded)",
    y = "Subcluster",
    x = "Gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted")
  )

ggsave(file.path(OUT_DIR, "01_DotPlot_Before_Annotation.pdf"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(OUT_DIR, "01_DotPlot_Before_Annotation.png"), p1, width = 8, height = 6, dpi = 300, bg = "white")

cat("\n Info...\n")

manual_map <- c(
  "0" = "Basal",
  "1" = "Club",
  "2" = "Luminal",
  "3" = "Basal",
  "4" = "Basal",
  "7" = "Luminal"
)

epi$epi_subtype <- NA_character_

for (s in names(manual_map)) {
  cells_in_sub <- colnames(epi)[epi$subcluster == s &
                                  epi$malignant_status %in% c("Normal-like", "Reference")]
  epi$epi_subtype[cells_in_sub] <- manual_map[s]
}

epi$epi_subtype[epi$subcluster == "5"] <- "Malignant_1"
epi$epi_subtype[epi$subcluster == "6"] <- "Malignant_2"
epi$epi_subtype[epi$malignant_status == "Malignant" & is.na(epi$epi_subtype)] <- "Malignant_Other"
epi$epi_subtype[is.na(epi$epi_subtype)] <- "Unknown"
epi$epi_subtype <- factor(epi$epi_subtype)

cat("Info (Info):\n")
print(table(epi$epi_subtype, useNA = "ifany"))

cat("\n Info Malignant_Other...\n")
keep_cells <- colnames(epi)[epi$epi_subtype != "Malignant_Other"]
epi <- subset(epi, cells = keep_cells)
epi$epi_subtype <- droplevels(epi$epi_subtype)

cat("Info:", ncol(epi), "cells\n")
print(table(epi$epi_subtype, useNA = "ifany"))

cat("\n Plotting Info (Info cells)...\n")

normal_anno_cells <- colnames(epi)[epi$epi_subtype %in% c("Basal", "Luminal", "Club")]
epi_normal_anno <- subset(epi, cells = normal_anno_cells)
epi_normal_anno$epi_subtype <- factor(epi_normal_anno$epi_subtype, levels = c("Basal", "Luminal", "Club"))
Idents(epi_normal_anno) <- "epi_subtype"

p2 <- DotPlot(epi_normal_anno, features = all_markers, cols = c("white", "#1f77b4"), dot.scale = 7) +
  coord_flip() +
  labs(
    title = "Marker Expression by Annotated Subtype",
    subtitle = "Normal epithelial cells only (Basal / Luminal / Club)",
    y = "Annotated Subtype",
    x = "Gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted")
  )

ggsave(file.path(OUT_DIR, "02_DotPlot_After_Annotation_NormalOnly.pdf"), p2, width = 7, height = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "02_DotPlot_After_Annotation_NormalOnly.png"), p2, width = 7, height = 5, dpi = 300, bg = "white")

cat("\n Plotting Info UMAP...\n")

if (!"umap" %in% names(epi@reductions)) {
  epi <- ScaleData(epi, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:20, verbose = FALSE)
}

epi$epi_subtype <- factor(epi$epi_subtype, levels = c("Basal", "Luminal", "Club", "Malignant_1", "Malignant_2", "Unknown"))

subtype_colors <- c(
  "Basal"       = "#1f77b4",
  "Luminal"     = "#2ca02c",
  "Club"        = "#9467bd",
  "Malignant_1" = "#E64B35",
  "Malignant_2" = "#FF6B6B",
  "Unknown"     = "#7f7f7f"
)

p3 <- DimPlot(epi, group.by = "epi_subtype", cols = subtype_colors,
              pt.size = 0.6, label = TRUE, repel = TRUE, label.size = 4) +
  ggtitle("Epithelial Subtype Annotation") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "right")

ggsave(file.path(OUT_DIR, "03_UMAP_After_Annotation.pdf"), p3, width = 10, height = 7, dpi = 300)
ggsave(file.path(OUT_DIR, "03_UMAP_After_Annotation.png"), p3, width = 10, height = 7, dpi = 300, bg = "white")

saveRDS(epi, file.path(OUT_DIR, "epi_all_annotated.rds"))

cat("\n All completed! Output:", OUT_DIR, "\n")
cat("01_DotPlot_Before_Annotation.pdf/png - Info (Info)\n")
cat("02_DotPlot_After_Annotation_NormalOnly.pdf/png - Info (Info)\n")
cat("03_UMAP_After_Annotation.pdf/png - Info UMAP(None Malignant_Other)\n")
cat("epi_all_annotated.rds - Info \n")
cat("========================================\n")

sessionInfo()
