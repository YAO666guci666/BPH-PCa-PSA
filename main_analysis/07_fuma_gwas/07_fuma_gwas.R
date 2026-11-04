if (!exists(".shared_pipeline_call")) rm(list = ls())

gc()

getwd()

.libPaths()

library(data.table)
library(Fast2TWAS)

if (!exists("inputFile")) inputFile <- c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")

for (trait in inputFile) {
  if (! dir.exists(trait)) {
    dir.create(trait)
  }
}

stable6 <- data.frame()
for (trait in inputFile) {

  GenomicRiskLoci_dir <- list.files(
    path = paste0("../07_fuma_gwas/",trait,"/"),
    pattern = "^FUMA.*",
    full.names = TRUE,
    include.dirs = TRUE
  ) %>% .[file.info(.)$isdir]

  GenomicRiskLoci=fread(paste0(GenomicRiskLoci_dir,"/GenomicRiskLoci.txt"))

  GenomicRiskLoci$chr <- gsub("^23$|^x$", "X", as.character(GenomicRiskLoci$chr))
  GenomicRiskLoci$Trait <- trait
  GenomicRiskLoci <- subset(GenomicRiskLoci, select=c(Trait,uniqID,rsID,pos))

  gwas_data <- fread(paste0("../02_harmonize_gwas/",trait,"_hg38_qc.txt"))

  gwas_data$CHR <- gsub("^23$|^x$", "X", as.character(gwas_data$CHR))

  pre_stable6 <- merge(GenomicRiskLoci, gwas_data, by.x="rsID", by.y="SNP")

  setnames(pre_stable6, "BP", "hg38")

  setnames(pre_stable6, "uniqID", "GWAS_locus")

  setnames(pre_stable6, "pos", "hg19")

  pre_stable6$CHR <- paste0("chr", pre_stable6$CHR)

  pre_stable6$GWAS_locus <- paste0("chr", pre_stable6$GWAS_locus)

  pre_stable6$GWAS_locus <- gsub("^chr23", "chrX", pre_stable6$GWAS_locus)

  stable6 <- rbind(stable6, pre_stable6)

  fwrite(stable6,
         file="stable6.txt",
         sep="\t")
}

sessionInfo()
