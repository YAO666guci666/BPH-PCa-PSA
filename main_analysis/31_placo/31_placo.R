library(QTLMR)

pairs <- list(
  c("BPH", "PSA"),
  c("BPH", "PCA"),
  c("PCA", "PSA")
)

for(pair in pairs){
  t1 <- pair[1]
  t2 <- pair[2]
  p_name <- paste0(t1, "_", t2)

  cat("\n>>> Running PLACO:", p_name, "\n")

  out_dir <- file.path("./PLACO_pairwise", p_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  dat_PLACO <- PLACO_analysis(
    GWASfile = c(paste0("../01_preprocess_gwas/", t1, ".txt"),
                 paste0("../01_preprocess_gwas/", t2, ".txt")),
    p_threshold = 5e-08,
    filter_Zscores = TRUE,
    FDR_method = "fdr",
    decor = FALSE,
    save_path = out_dir,
    cores = 4
  )

  cat("Completed:", p_name, "\n")
  cat("Output path:", out_dir, "\n")
}

cat("\n========================================\n")
cat("All three PLACO All completed \n")
cat("========================================\n")

library(ggplot2)
library(dplyr)

pair_map <- data.frame(
  label = c("BPH_PSA", "BPH_PCa", "PSA_PCa"),
  dir   = c("BPH_PSA", "BPH_PCA", "PCA_PSA")
)

all_placo <- lapply(1:nrow(pair_map), function(i) {
  p_label <- pair_map$label[i]
  p_dir   <- pair_map$dir[i]

  f <- file.path("./PLACO_pairwise", p_dir, "placo_dat.txt")
  cat("Reading:", p_label, "->", f, "\n")

  if(!file.exists(f)) stop(paste("File not found:", f))

  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  df$pair <- p_label
  return(df)
}) %>% bind_rows()

colnames(all_placo) <- tolower(colnames(all_placo))
if("p" %in% colnames(all_placo) & !"p_placo" %in% colnames(all_placo)) all_placo$p_placo <- all_placo$p
if("p.placo" %in% colnames(all_placo) & !"p_placo" %in% colnames(all_placo)) all_placo$p_placo <- all_placo$p.placo
if("placo_p" %in% colnames(all_placo) & !"p_placo" %in% colnames(all_placo)) all_placo$p_placo <- all_placo$placo_p

chr_lengths <- all_placo %>%
  group_by(pair, chr) %>%
  summarise(max_pos = max(as.numeric(pos), na.rm = TRUE), .groups = "drop") %>%
  group_by(pair) %>%
  mutate(offset = c(0, cumsum(head(max_pos, -1) + 1e8))) %>%
  ungroup()

plot_data <- all_placo %>%
  left_join(chr_lengths, by = c("pair", "chr")) %>%
  mutate(
    pos       = as.numeric(pos) + offset,
    logp      = -log10(as.numeric(p_placo)),
    is_sig    = as.numeric(p_placo) < 5e-8,
    chr       = as.integer(chr),
    chr_color = ifelse(chr %% 2 == 1, "odd", "even")
  ) %>%
  filter(!is.na(pos), !is.na(logp), chr %in% 1:22)

axis_df <- plot_data %>%
  group_by(pair, chr) %>%
  summarise(center = mean(pos), .groups = "drop") %>%
  filter(chr <= 22)

p <- ggplot(plot_data, aes(x = pos, y = logp)) +
  geom_point(data = subset(plot_data, !is_sig),
             aes(color = chr_color, size = chr_color),
             alpha = 0.6) +
  geom_point(data = subset(plot_data, is_sig),
             color = "#E74C3C", size = 1.5, alpha = 0.9) +
  geom_hline(yintercept = -log10(5e-8),
             color = "#E74C3C", linetype = "dashed", linewidth = 0.3) +
  facet_grid(pair ~ ., scales = "free_y", space = "free_x") +
  scale_color_manual(values = c("odd" = "#9B8EC2", "even" = "#7EC8B8"), guide = "none") +
  scale_size_manual(values = c("odd" = 1.2, "even" = 0.8), guide = "none") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr,
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0, 0.5)) +
  labs(x = "Chromosome", y = expression(-log[10](P)),
       title = "PLACO Multi-trait Manhattan Plot") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

out_png <- "./PLACO_pairwise/PLACO_Manhattan_3pairs.png"
ggsave(out_png, p, width = 12, height = 10, dpi = 300, bg = "white")
cat("Figure saved:", out_png, "\n")

library(ggplot2)
library(dplyr)

pair_map <- data.frame(
  label = c("BPH_PSA", "BPH_PCa", "PSA_PCa"),
  dir   = c("BPH_PSA", "BPH_PCA", "PCA_PSA")
)

all_placo <- lapply(1:nrow(pair_map), function(i) {
  p_label <- pair_map$label[i]
  p_dir   <- pair_map$dir[i]
  f <- file.path("./PLACO_pairwise", p_dir, "placo_dat.txt")
  cat("Reading placo_dat:", p_label, "\n")
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  df$pair <- p_label
  return(df)
}) %>% bind_rows()

colnames(all_placo) <- tolower(colnames(all_placo))
if("p.placo" %in% colnames(all_placo) & !"p_placo" %in% colnames(all_placo)) all_placo$p_placo <- all_placo$p.placo

read_lead <- function(p_label, p_dir) {

  f <- file.path("./PLACO_pairwise", p_dir, p_dir, "leadSNPs.txt")

  if(!file.exists(f)) {
    cat("Not found leadSNPs.txt, trying IndSigSNPs.txt:", p_label, "\n")
    f <- file.path("./PLACO_pairwise", p_dir, p_dir, "IndSigSNPs.txt")
  }

  if(!file.exists(f)) {
    cat("Warning: No independent-locus file for", p_label, "\n")
    return(data.frame(snp = character(0), pair = character(0), stringsAsFactors = FALSE))
  }

  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- tolower(colnames(df))

  id_col <- intersect(c("rsid", "snp", "rs", "markername"), colnames(df))
  if(length(id_col) == 0) id_col <- colnames(df)[1]

  res <- data.frame(snp = as.character(df[[id_col[1]]]), pair = p_label, stringsAsFactors = FALSE)
  cat("Read", nrow(res), "lead SNPs for", p_label, "\n")
  return(res)
}

lead_snps <- lapply(1:nrow(pair_map), function(i) {
  read_lead(pair_map$label[i], pair_map$dir[i])
}) %>% bind_rows()

all_placo$is_lead <- paste(all_placo$snp, all_placo$pair) %in% paste(lead_snps$snp, lead_snps$pair)

chr_lengths <- all_placo %>%
  group_by(pair, chr) %>%
  summarise(max_pos = max(as.numeric(pos), na.rm = TRUE), .groups = "drop") %>%
  group_by(pair) %>%
  mutate(offset = c(0, cumsum(head(max_pos, -1) + 1e8))) %>%
  ungroup()

plot_data <- all_placo %>%
  left_join(chr_lengths, by = c("pair", "chr")) %>%
  mutate(
    pos       = as.numeric(pos) + offset,
    logp      = -log10(as.numeric(p_placo)),
    chr       = as.integer(chr),
    chr_color = ifelse(chr %% 2 == 1, "odd", "even")
  ) %>%
  filter(!is.na(pos), !is.na(logp), chr %in% 1:22)

axis_df <- plot_data %>%
  group_by(pair, chr) %>%
  summarise(center = mean(pos), .groups = "drop") %>%
  filter(chr <= 22)

p <- ggplot(plot_data, aes(x = pos, y = logp)) +

  geom_point(data = subset(plot_data, !is_lead),
             aes(color = chr_color, size = chr_color),
             alpha = 0.6) +

  geom_point(data = subset(plot_data, is_lead),
             color = "#E74C3C", size = 2.2, alpha = 0.95) +

  geom_hline(yintercept = -log10(5e-8),
             color = "#E74C3C", linetype = "dashed", linewidth = 0.3) +
  facet_grid(pair ~ ., scales = "free_y", space = "free_x") +
  scale_color_manual(values = c("odd" = "#9B8EC2", "even" = "#7EC8B8"), guide = "none") +
  scale_size_manual(values = c("odd" = 1.2, "even" = 0.8), guide = "none") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr,
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0, 0.5)) +
  labs(x = "Chromosome", y = expression(-log[10](P)),
       title = "PLACO Multi-trait Manhattan Plot (Lead SNPs Highlighted)") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

out_png <- "./PLACO_pairwise/PLACO_Manhattan_3pairs_lead.png"
ggsave(out_png, p, width = 12, height = 10, dpi = 300, bg = "white")
cat("Figure saved:", out_png, "\n")

sessionInfo()
