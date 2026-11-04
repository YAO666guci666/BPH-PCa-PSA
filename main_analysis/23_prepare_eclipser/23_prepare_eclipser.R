rm(list=ls())

gc()

getwd()

.libPaths()

trait="PC"
df <- read.csv(paste0(trait, ".csv"), header = TRUE, stringsAsFactors = FALSE)

colnames(df)[colnames(df) == "cell_type"] <- "celltype"
colnames(df)[colnames(df) == "p_val_adj"] <- "pvals_fdr"
colnames(df)[colnames(df) == "avg_log2FC"] <- "log2FC"

new_df <- df[, c("celltype", "gene", "tissue", "pvals_fdr", "log2FC")]

head(new_df)

output_filename <- paste0(trait, "_train.csv")
write.csv(new_df, output_filename, row.names = FALSE)

sessionInfo()
