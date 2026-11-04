library(dplyr)

convert_sym <- function(infile, outfile) {
  dat <- read.table(infile, header = TRUE, sep = "", comment.char = "#", stringsAsFactors = FALSE)
  stopifnot(all(c("SYMBOL", "ZSTAT") %in% colnames(dat)))
  res <- dat %>%
    select(GENE = SYMBOL, ZSTAT) %>%
    mutate(ZSTAT = as.numeric(ZSTAT)) %>%
    filter(!is.na(ZSTAT)) %>%
    group_by(GENE) %>%
    slice_max(order_by = abs(ZSTAT), n = 1, with_ties = FALSE) %>%
    ungroup()
  colnames(res) <- c("GENE", "zscore")
  write.table(res, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
  cat(outfile, "Completed, Total", nrow(res), "genes \n")
}

convert_6319 <- function(csv_file, outfile) {
  stopifnot(file.exists(csv_file))
  dat <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE,
                  check.names = FALSE, fileEncoding = "UTF-8-BOM")
  stopifnot(all(c("Gene", "ZSTAT") %in% colnames(dat)))
  res <- dat %>%
    select(GENE = Gene, ZSTAT) %>%
    mutate(ZSTAT = as.numeric(ZSTAT)) %>%
    filter(!is.na(ZSTAT)) %>%
    group_by(GENE) %>%
    slice_max(order_by = abs(ZSTAT), n = 1, with_ties = FALSE) %>%
    ungroup()
  colnames(res) <- c("GENE", "zscore")
  write.table(res, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
  cat(outfile, "Completed, Total", nrow(res), "genes \n")
}

convert_sym("PCa_PSA_magma.genes.out", "PCa_PSA.tsv")
convert_sym("BPH_PSA_magma.genes.out", "BPH_PSA.tsv")

convert_6319("PCA_PSA_6319_MAGMA_all_genes.csv", "PCA_PSA_6319.tsv")
convert_6319("BPH_PSA_6319_MAGMA_all_genes.csv", "BPH_PSA_6319.tsv")
