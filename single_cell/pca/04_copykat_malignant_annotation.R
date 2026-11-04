rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 100)
set.seed(1234)

library(Seurat)
library(copykat)
library(tidyverse)
library(ggplot2)

if (!requireNamespace("conflicted", quietly = TRUE)) install.packages("conflicted")
library(conflicted)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)

EPI_PATH <- "01_Epithelial_Subclustering/13_epithelial_subclustered.rds"

NORMAL_GROUP_VALUE <- "Healthy"

TUMOR_GROUPS <- c("Tumor_LowRisk", "Tumor_HighRisk", "Tumor")

OUT_DIR <- "02_CopyKAT_Combined_Tumor"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

NCORES <- 8
NGENE_CHR <- 5

cat("Info data...\n")
epi_all <- readRDS(EPI_PATH)

cat("raw Info group Info:\n")
print(table(epi_all$group, useNA = "ifany"))

cat("\n Check Merging Info group...\n")

epi_all$group <- as.character(epi_all$group)

epi_all$major_group <- ifelse(
  epi_all$group %in% TUMOR_GROUPS,
  "Tumor",
  epi_all$group
)

epi_all$major_group <- factor(epi_all$major_group, levels = c("Healthy", "Tumor"))

cat("Merging Info group Info:\n")
print(table(epi_all$major_group, useNA = "ifany"))

epi_all$group <- epi_all$major_group

normal_cells <- colnames(epi_all)[epi_all$major_group == NORMAL_GROUP_VALUE]
cat("\n Info cells (Healthy):", length(normal_cells), "\n")

if (length(normal_cells) == 0) {
  stop("Not found Info cells, check NORMAL_GROUP_VALUE Info")
}

cat("\n Preparing Info...\n")
DefaultAssay(epi_all) <- "RNA"
exp_raw <- as.matrix(GetAssayData(epi_all, layer = "counts"))
exp_raw <- exp_raw[rowSums(exp_raw) > 0, ]

cat("Info:", nrow(exp_raw), "genes x", ncol(exp_raw), "cells \n")

cat("\n Info rows CopyKAT(Info group)...\n")
cat("Note: Info 10~30 Info...\n")

copykat_result <- copykat(
  rawmat = exp_raw,
  id.type = "S",
  cell.line = "no",
  ngene.chr = NGENE_CHR,
  win.size = 25,
  KS.cut = 0.1,
  norm.cell.names = normal_cells,
  n.cores = NCORES,
  output.seg = FALSE,
  plot.genes = FALSE
)

saveRDS(copykat_result, file.path(OUT_DIR, "copykat_result.rds"))
cat("CopyKAT raw results Info Saved \n")

cat("\n Info CopyKAT results...\n")

pred_df <- data.frame(
  cell_id = copykat_result$prediction$cell.names,
  copykat_raw = copykat_result$prediction$copykat.pred,
  stringsAsFactors = FALSE
)

pred_df$copykat_class <- ifelse(
  pred_df$copykat_raw %in% c("diploid", "Diploid"), "Diploid",
  ifelse(pred_df$copykat_raw %in% c("aneuploid", "Aneuploid"), "Aneuploid", NA)
)

rownames(pred_df) <- pred_df$cell_id
epi_all$copykat_class <- pred_df[colnames(epi_all), "copykat_class"]

cat("CopyKAT Info:\n")
print(table(epi_all$copykat_class, useNA = "ifany"))

cat("\n Info CopyKAT Info rows Info...\n")

epi_all$malignant_status <- dplyr::case_when(
  epi_all$major_group == NORMAL_GROUP_VALUE ~ "Reference",
  epi_all$copykat_class == "Aneuploid" ~ "Malignant",
  epi_all$copykat_class == "Diploid" ~ "Normal-like",
  TRUE ~ NA_character_
)

cat("Info:\n")
print(table(epi_all$malignant_status, epi_all$major_group, useNA = "ifany"))

saveRDS(epi_all, file.path(OUT_DIR, "epi_all_copykat_annotated.rds"))

cat("\n Plotting UMAP...\n")

if (!"umap" %in% names(epi_all@reductions)) {
  cat("Info rows UMAP...\n")
  epi_all <- NormalizeData(epi_all) %>% FindVariableFeatures() %>% ScaleData() %>%
    RunPCA() %>% RunUMAP(dims = 1:20)
}

copykat_colors <- c("Diploid" = "#4DBBD5", "Aneuploid" = "#E64B35")

p_umap_ck <- DimPlot(epi_all, group.by = "copykat_class", cols = copykat_colors,
                     pt.size = 0.5, order = c("Aneuploid", "Diploid")) +
  labs(title = "CopyKAT (Combined Tumor vs Healthy Ref)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "01_UMAP_CopyKAT.pdf"), p_umap_ck, width = 8, height = 7)
ggsave(file.path(OUT_DIR, "01_UMAP_CopyKAT.png"), p_umap_ck, width = 8, height = 7, dpi = 300)

gold_colors <- c("Reference" = "#4DBBD5", "Normal-like" = "#91D1C2", "Malignant" = "#E64B35")

p_umap_gold <- DimPlot(epi_all, group.by = "malignant_status", cols = gold_colors,
                       pt.size = 0.5, order = c("Malignant", "Normal-like", "Reference")) +
  labs(title = "Gold Standard: Malignant vs Normal-like") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "02_UMAP_GoldStandard.pdf"), p_umap_gold, width = 8, height = 7)

cat("Plotting Info group Info...\n")

group_stats <- epi_all@meta.data %>%
  dplyr::filter(!is.na(copykat_class)) %>%
  dplyr::count(major_group, copykat_class, name = "n") %>%
  dplyr::group_by(major_group) %>%
  dplyr::mutate(pct = n / sum(n) * 100)

p_stack_group <- ggplot(group_stats, aes(x = major_group, y = pct, fill = copykat_class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = copykat_colors) +
  labs(
    title = "CopyKAT Classification: Healthy vs Combined Tumor",
    x = "Group",
    y = "Percentage (%)",
    fill = "CopyKAT"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "03_StackedBar_Group.pdf"), p_stack_group, width = 6, height = 6)

cat("Plotting Info (Merging Info group)...\n")

sub_stats <- epi_all@meta.data %>%
  dplyr::filter(
    major_group == "Tumor",
    !is.na(malignant_status),
    !is.na(subcluster)
  ) %>%
  dplyr::count(subcluster, malignant_status, name = "n") %>%
  dplyr::group_by(subcluster) %>%
  dplyr::mutate(pct = n / sum(n) * 100)

sub_stats$subcluster <- factor(
  sub_stats$subcluster,
  levels = sort(as.numeric(as.character(unique(sub_stats$subcluster))))
)

p_stack_sub <- ggplot(sub_stats, aes(x = subcluster, y = pct, fill = malignant_status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = ifelse(pct > 5, paste0(round(pct, 1), "%"), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = gold_colors) +
  labs(
    title = "Malignant Proportion by Subcluster (Combined Tumor)",
    x = "Subcluster",
    y = "Percentage (%)",
    fill = "Status"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(OUT_DIR, "04_StackedBar_Subcluster_Tumor.pdf"),
       p_stack_sub, width = 8, height = 6)

cat("\n Output Info...\n")

stats_summary <- epi_all@meta.data %>%
  dplyr::filter(!is.na(malignant_status)) %>%
  dplyr::count(major_group, subcluster, malignant_status, name = "n") %>%
  tidyr::pivot_wider(names_from = malignant_status, values_from = n, values_fill = 0)

status_cols <- setdiff(names(stats_summary), c("major_group", "subcluster"))
stats_summary$total <- rowSums(stats_summary[, status_cols, drop = FALSE], na.rm = TRUE)

if ("Malignant" %in% status_cols) {
  stats_summary$malignant_pct <- round(stats_summary$Malignant / stats_summary$total * 100, 2)
} else {
  stats_summary$malignant_pct <- 0
}

stats_summary <- stats_summary %>% dplyr::arrange(major_group, desc(malignant_pct))
write.csv(stats_summary, file.path(OUT_DIR, "GoldStandard_Stats.csv"), row.names = FALSE)

cat("Info Saved:", file.path(OUT_DIR, "GoldStandard_Stats.csv"), "\n")
cat("\n========================================\n")
cat("CopyKAT Info Completed!\n")
cat("========================================\n")

sessionInfo()
