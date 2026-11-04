rm(list=ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)
library(data.table)
library(rtracklayer)

inputFile=c("ieu-b-4809_PC","GCST90476277_PSA","finngen_BPH")

ref_pop="g1000_eur"

for (trait in inputFile) {

  stable6 <- fread(paste0("../07_fuma_gwas/stable6.txt"))

  stable6_filter <- subset(stable6, Trait == trait)

  gwas_file <- paste0("../02_harmonize_gwas/",trait,"_hg38_qc.txt")
  gwas_data <- fread(gwas_file)

  GenomicRiskLoci_dir <- list.files(
    path = paste0("../07_fuma_gwas/",trait,"/"),
    pattern = "^FUMA.*",
    full.names = TRUE,
    include.dirs = TRUE
  ) %>% .[file.info(.)$isdir]

  GenomicRiskLoci=fread(paste0(GenomicRiskLoci_dir,"/GenomicRiskLoci.txt"))

  if (!"eaf" %in% colnames(gwas_data)) {
    new_gwas_file <- paste0("../rawdata/sumstats_rehab/result/",trait,".txt.rehabed.tsv")
    gwas_data_filter <- fread(new_gwas_file)[, c("SNP", "effect_allele", "other_allele","EAF_rehab", "beta", "se", "pval","N")]
  }else{

  gwas_data_filter <- fread(gwas_file)[, c("SNP", "effect_allele", "other_allele","eaf", "beta", "se", "pval","N")]
  }

  names(gwas_data_filter) <- c("SNP", "A1", "A2", "freq", "b", "se", "p","N")

  gwas_data_filter <- gwas_data_filter[gwas_data_filter$freq != ".", ]

  fwrite(gwas_data_filter,
         file = paste0(trait,".ma"),
         sep = "\t")

  for (snp in stable6_filter$rsID) {

  trait_snp_dir <- paste0("./result/",trait,"/",snp,"/")
  if (!dir.exists(trait_snp_dir)) {
    dir.create(trait_snp_dir, recursive  = TRUE)
  }
  chr <- gsub("chr", "", subset(stable6_filter, rsID == snp)$CHR)
  snplist_file <- paste0(trait_snp_dir,snp,"_cond.snplist")
  write(x = snp, file = snplist_file, append = F)

  GenomicRiskLoci_snp <- subset(GenomicRiskLoci,rsID==snp)

  gr_hg19 <- GRanges(
    seqnames = paste0("chr", GenomicRiskLoci_snp$chr),
    ranges = IRanges(start = c(GenomicRiskLoci_snp$start,GenomicRiskLoci_snp$end), width = 1)
  )

  chain <- import.chain("../reference/hg19ToHg38.over.chain")
  gr_hg38 <- unlist(liftOver(gr_hg19, chain))

  gwas_data_hg38 <- data.frame(
    CHR = gsub("chr","",seqnames(gr_hg38)),
    BP_hg38 = start(gr_hg38)
  )

  CHR <- unique(gwas_data_hg38$CHR)
  start <- min(gwas_data_hg38$BP_hg38)
  end <- max(gwas_data_hg38$BP_hg38)

  gwas_region <- gwas_data %>% subset(.,CHR == CHR & BP >= start & BP <= end)

  snp_list <- gwas_region$SNP

  extract_file <- paste0(snplist_file,trait,"_",snp,"_extract.txt")

  write(x = snp_list, file = extract_file, append = F)

  gcta_path <- friendlypostGWASsoft::gcta_binary()
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd","sh")

  GCTA_cojo_cond <- paste0(shQuote(gcta_path, type = shell),
                  " --bfile ","../reference/magma/g1000/",ref_pop,
                  " --chr ",CHR,
                  " --extract ",extract_file,
                  " --maf 0.001 ",
                  " --cojo-file ",paste0(trait,".ma"),
                  " --cojo-cond ",snplist_file,
                  " --out ",paste0(trait_snp_dir,trait,"_cond_",snp))

  system(GCTA_cojo_cond)

  }
}

sessionInfo()
