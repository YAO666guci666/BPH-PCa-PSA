rm(list=ls())
gc()

library(Fast2TWAS)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(stringr)
library(patchwork)

inputFile <- c("ieu-b-4809_PC", "finngen_BPH", "GCST90476277_PSA")
tissue_target <- "Prostate"

combined_data <- data.frame()

for (trait in inputFile) {
  enrich_file <- paste0("../15_process_gene_enrich/stable_", trait, ".GSEA_QTLEnrich.txt")

  if (file.exists(enrich_file)) {
    enrich_data <- fread(enrich_file)
    enrich_data_subset <- subset(enrich_data, Tissue == tissue_target & QTL_type %in% c("eQTL", "sQTL"))

    if (nrow(enrich_data_subset) > 0) {
      enrich_data_subset$Trait <- trait
      combined_data <- bind_rows(combined_data, enrich_data_subset)
    }
  } else {
    warning("File not found:", enrich_file)
  }
}

if (nrow(combined_data) == 0) {
  stop("combined_data is empty, check the file paths and names.")
}

combined_df <- as.data.frame(combined_data)

combined_df <- combined_df[combined_df$resource != "MGI", , drop = FALSE]

combined_df$FDR_db <- NA_real_

for (tr in unique(combined_df$Trait)) {
  for (res in unique(combined_df$resource)) {
    for (qt in unique(combined_df$QTL_type)) {
      idx <- combined_df$Trait == tr & combined_df$resource == res & combined_df$QTL_type == qt
      if (sum(idx) > 0) {
        combined_df$FDR_db[idx] <- p.adjust(combined_df$empirical_pval[idx], method = "BH")
      }
    }
  }
}

m_df <- aggregate(empirical_pval ~ Trait + resource + QTL_type,
                  data = combined_df, FUN = length)
names(m_df)[names(m_df) == "empirical_pval"] <- "m"

k_df <- aggregate(FDR_db < 0.05 ~ Trait + resource + QTL_type,
                  data = combined_df, FUN = sum)
names(k_df)[names(k_df) == "FDR_db < 0.05"] <- "k"

line_df <- merge(m_df, k_df, by = c("Trait", "resource", "QTL_type"))
line_df$x_line <- -log10(0.05 * pmax(line_df$k, 1) / line_df$m)

sig_hits <- combined_df[combined_df$FDR_db < 0.05,
                        c("Trait", "resource", "QTL_type", "gene_set", "empirical_pval", "FDR_db")]
cat("Database FDR<0.05 pathways:\n")
if (nrow(sig_hits) > 0) {
  print(sig_hits[order(sig_hits$FDR_db), ])
} else {
  cat("(None)\n")
}
cat("Tests per stratum m and passed k:\n")
print(line_df)

combined_data <- combined_df

plot_data <- as.data.frame(combined_data[, c("Trait", "resource", "QTL_type", "gene_set", "empirical_pval")])
names(plot_data) <- c("Trait", "exp", "QTL_type", "Pathway_simplified", "P")

plot_data_top <- plot_data %>%
  group_by(Trait, exp, QTL_type) %>%
  slice_min(P, n = 3, with_ties = FALSE) %>%
  ungroup()

dd <- as.data.frame(plot_data_top)
dd$Pathway_simplified <- gsub("_", " ", tolower(dd$Pathway_simplified))

dd$Pathway_simplified <- str_wrap(dd$Pathway_simplified, width = 55)

dd$Trait <- factor(dd$Trait, levels = inputFile)
dd$QTL_type <- factor(dd$QTL_type, levels = c("eQTL", "sQTL"))
dd$exp <- factor(dd$exp)
dd$logP <- (-1) * log10(dd$P)

if (!dir.exists("results")) {
  dir.create("results")
}

trait_mapping <- c(
  "ieu-b-4809_PC" = "PCa",
  "finngen_BPH" = "BPH",
  "GCST90476277_PSA" = "PSA"
)

exp_block_A <- c("GO_BP", "GO_CC", "GO_MF")
exp_block_B <- c("HALLMARK", "KEGG", "REACTOME")

for (i in 1:length(inputFile)) {
  trait_name <- inputFile[i]
  trait_simple <- trait_mapping[trait_name]

  dd_trait <- dd[dd$Trait == trait_name, , drop = FALSE]
  if (nrow(dd_trait) == 0) next

  line_trait <- line_df[line_df$Trait == trait_name, , drop = FALSE]
  line_trait$exp <- line_trait$resource

  block_plots <- list()
  block_widths <- c()
  block_heights <- c()

  for (block_id in c("A", "B")) {
    if (block_id == "A") {
      exps_use <- exp_block_A
    } else {
      exps_use <- exp_block_B
    }

    exps_use <- intersect(exps_use, unique(as.character(dd_trait$exp)))
    if (length(exps_use) == 0) next

    dd_block <- dd_trait[as.character(dd_trait$exp) %in% exps_use, , drop = FALSE]
    dd_block$exp <- factor(as.character(dd_block$exp), levels = exps_use)

    line_block <- line_trait[as.character(line_trait$exp) %in% exps_use, , drop = FALSE]

    max_pathway_length <- max(vapply(
      strsplit(as.character(dd_block$Pathway_simplified), "\n"),
      function(l) max(nchar(l)), numeric(1)
    ))

    n_exp_block <- length(exps_use)
    block_height <- nrow(dd_block) * 2.0 + n_exp_block * 1.6 + 1.6

    block_width <- 15 + max_pathway_length * 0.46

    if (block_id == "A") {
      margin_use <- unit(c(3, 0.3, 3, 6), "mm")
    } else {
      margin_use <- unit(c(3, 3, 3, 6), "mm")
    }

    if (nrow(line_block) > 0) {
      x_max <- max(c(max(dd_block$logP) * 1.08, line_block$x_line * 1.03))
    } else {
      x_max <- max(dd_block$logP) * 1.08
    }

    p_block <- ggplot(data = dd_block, aes(y = fct_reorder(Pathway_simplified, logP), x = logP)) +
      geom_bar(stat = "identity", fill = "#7C9895", width = 0.65, alpha = 0.85) +
      labs(y = "", x = expression("-log"[10]*"(P value)"), color = NULL, linetype = NULL) +

    facet_grid(exp ~ QTL_type, scales = "free_y", space = "free_y") +

    geom_vline(aes(xintercept = -log10(0.05),
                   color = "Nominal P = 0.05",
                   linetype = "Nominal P = 0.05"),
               alpha = 0.7, linewidth = 1.0) +

    geom_vline(data = line_block,
               aes(xintercept = x_line,
                   color = "Within-database FDR = 0.05",
                   linetype = "Within-database FDR = 0.05"),
               alpha = 0.7, linewidth = 1.0) +
      scale_color_manual(values = c("Nominal P = 0.05" = "red",
                                    "Within-database FDR = 0.05" = "blue")) +
      scale_linetype_manual(values = c("Nominal P = 0.05" = "dashed",
                                       "Within-database FDR = 0.05" = "dashed")) +

    guides(color = guide_legend(nrow = 1, order = 1),
           linetype = guide_legend(nrow = 1, order = 1)) +
      theme_bw() +
      theme(
        legend.position = "top",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "center",
        legend.spacing.x = unit(0.6, "cm"),
        legend.text = element_text(size = 21, face = "bold"),
        legend.key.width = unit(1.8, "cm"),
        legend.margin = margin(b = 0),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_line(colour = "gray85", linewidth = 0.15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),

        axis.text.y = element_text(size = 27, face = "bold", hjust = 1,
                                   lineheight = 1.0, margin = margin(r = 4)),
        axis.text.x = element_text(size = 22, face = "bold"),
        axis.title = element_text(size = 27, face = "bold"),
        strip.text = element_text(size = 24, face = "bold", margin = margin(1.8, 1.8, 1.8, 1.8, "mm")),
        plot.title = element_text(hjust = 0.5, size = 33, face = "bold", margin = margin(b = 3)),
        plot.margin = margin_use,
        panel.spacing = unit(3.5, "mm"),
        strip.background = element_rect(fill = "gray95", linewidth = 0.1)
      ) +
      coord_cartesian(xlim = c(0, x_max))

    block_plots[[block_id]] <- p_block
    block_widths[block_id] <- block_width
    block_heights[block_id] <- block_height
  }

  if (length(block_plots) == 0) next

  save_height <- max(block_heights)
  save_width <- sum(block_widths) + 2.5

  if (length(block_plots) == 2) {
    p <- (block_plots[["A"]] | block_plots[["B"]]) + plot_layout(widths = unname(block_widths[c("A", "B")]))
  } else {
    p <- block_plots[[1]]
  }
  p <- p + plot_annotation(
    title = trait_simple,
    theme = theme(plot.title = element_text(hjust = 0.5, size = 33, face = "bold"))
  )

  ggsave(
    paste0("./results/", trait_simple, "_CompletePathwayNames.png"),
    p,
    width = save_width,
    height = save_height,
    units = "cm",
    dpi = 300,
    limitsize = FALSE
  )

  cat("Saved", trait_simple, "figure, size:", round(save_width, 1), "x", round(save_height, 1), "cm\n")
}

cat("All analyses completed!\n")

sessionInfo()
