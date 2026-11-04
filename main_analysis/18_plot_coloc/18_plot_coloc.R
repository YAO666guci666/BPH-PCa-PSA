rm(list=ls())

gc()

getwd()

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(patchwork)

inputFile <- c("finngen_BPH", "GCST90476277_PSA", "ieu-b-4809_PC")
trait_labels <- c("finngen_BPH" = "BPH",
                  "GCST90476277_PSA" = "PSA",
                  "ieu-b-4809_PC" = "PCa")

cat_colors <- c("eGenes" = "#1B3E76",
                "sGenes" = "#FED71A",
                "Union"  = "#7C9895")

ref_theme <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18,
                              margin = margin(b = 12)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.5),
    axis.line.x = element_line(color = "grey60", linewidth = 0.5),
    axis.text = element_text(size = 14, color = "grey20"),
    axis.title.y = element_text(size = 14, color = "grey20",
                                margin = margin(r = 8)),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.7, "cm"),
    plot.margin = margin(15, 20, 12, 12)
  )

data <- fread("../17_prepare_coloc_plot/all_Num_GWAS_loci_coloc.txt")

data <- data %>% filter(Trait %in% inputFile)
data <- data[!data$Trait %in% c("Average", "Total loci"), ]

plot_left <- data.frame(
  Trait = rep(data$Trait, each = 3),
  Category = rep(c("eGenes", "sGenes", "Union"), times = nrow(data)),
  Count = c(rbind(data$loci_greater_than_1_colocalizing_eGene_union,
                  data$loci_greater_than_1_colocalizing_sGene_union,
                  data$loci_greater_than_1_colocalizing_eGene_sGene_union))
)
plot_left$Trait <- factor(trait_labels[plot_left$Trait], levels = c("BPH", "PSA", "PCa"))
plot_left$Category <- factor(plot_left$Category, levels = c("eGenes", "sGenes", "Union"))

print(plot_left)

p1 <- ggplot(plot_left, aes(x = Trait, y = Count, fill = Category)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.78),
            vjust = -0.6, size = 5.2, fontface = "bold", color = "grey15") +
  scale_fill_manual(values = cat_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "GWAS loci with >=1 colocalizing e/sQTLs",
       x = NULL, y = "Number of loci") +
  ref_theme +
  theme(legend.position = "top")

dat <- fread("../17_prepare_coloc_plot/all_NumQTL_coloc_locus.txt")

dat <- dat %>% filter(eGene_sGene_union > 0)

dat <- dat %>% filter(Trait %in% inputFile)

summary_stats <- dat %>%
  pivot_longer(
    cols = c(eGenes, sGenes, eGene_sGene_union),
    names_to = "Category",
    values_to = "Count"
  ) %>%
  mutate(Category = recode(Category,
                           "eGenes" = "eGenes",
                           "sGenes" = "sGenes",
                           "eGene_sGene_union" = "Union")) %>%
  group_by(Trait, Category) %>%
  summarize(
    Mean = mean(Count),
    SD = sd(Count),
    .groups = "drop"
  )
summary_stats$Trait <- factor(trait_labels[summary_stats$Trait], levels = c("BPH", "PSA", "PCa"))
summary_stats$Category <- factor(summary_stats$Category, levels = c("eGenes", "sGenes", "Union"))

summary_stats$lab_h <- ifelse(summary_stats$Category == "eGenes", 1.4,
                              ifelse(summary_stats$Category == "Union", -0.4, 0.5))
summary_stats$lab_v <- ifelse(summary_stats$Category == "sGenes", -1.1, 0.5)

print(summary_stats)

dodge <- position_dodge(width = 0.6)
p2 <- ggplot(summary_stats, aes(x = Trait, y = Mean, color = Category)) +
  geom_errorbar(aes(ymin = pmax(Mean - SD, 0), ymax = Mean + SD),
                width = 0.15, linewidth = 0.7, alpha = 0.85, position = dodge) +
  geom_point(size = 4, position = dodge) +
  geom_text(aes(label = sprintf("%.2f", Mean), hjust = lab_h, vjust = lab_v),
            position = dodge, size = 4.6,
            fontface = "bold", color = "grey25", show.legend = FALSE) +
  scale_color_manual(values = cat_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Colocalizing genes per locus",
       x = NULL, y = "Genes per locus (mean +/- SD)") +
  ref_theme +
  theme(legend.position = "none")

p <- p1 + p2 + plot_layout(ncol = 2)

print(p)

ggsave("Fig_2A_coloc_overview.pdf", p, width = 15, height = 6.5, dpi = 300)

sessionInfo()
