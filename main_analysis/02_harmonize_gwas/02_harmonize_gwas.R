if (!exists(".shared_pipeline_call")) rm(list = ls())
gc()
library(Fast2TWAS)
library(data.table)
if (!exists("inputFile")) inputFile <- c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")
maf=0.01
build="hg38"
for (trait in inputFile) {
  if (build=="hg19") {
    data=gwas_liftover(gwasfile=paste0("../01_preprocess_gwas/",trait,".txt"),
                       build_from="hg19",build_to="hg38",savefile=paste0(trait,"_hg38.txt"))
  } else if(build=="hg38") {
    data=fread(paste0("../01_preprocess_gwas/",trait,".txt"))
  }
  stopifnot(all(c("CHR","BP","SNP","effect_allele","other_allele") %in% names(data)))
  cat(trait,": before QC",nrow(data),"rows")
  data=subset(data,!(CHR==6 & BP > 25000000 & BP < 35000000))
  if("eaf" %in% names(data)) data=subset(data,eaf > maf & eaf < 1-maf)
  cat(trait,": after MHC and MAF filtering",nrow(data),"rows")
  fuma_data=subset(data,select=-c(CHR,BP))
  fuma_data=subset(fuma_data,!is.na(SNP) & SNP != "" & SNP != ".")
  fwrite(fuma_data,file=paste0(trait,"_hg38_qc_fuma.txt"),sep="	",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")
  system(paste("gzip -cf",paste0(trait,"_hg38_qc_fuma.txt"),">",paste0(trait,"_hg38_qc_fuma.gz")))
  fwrite(data,file=paste0(trait,"_hg38_qc.txt"),sep="	",quote=FALSE,row.names=FALSE)
  qtl_data=copy(data)
  setnames(qtl_data,c("pval","CHR","BP","other_allele"),c("gwas_p_value","chr","pos","non_effect_allele"),skip_absent=TRUE)
  qtl_data$chr=paste0("chr",qtl_data$chr)
  fwrite(qtl_data,file=paste0(trait,"_hg38_qc_QTLEnrich.txt"),sep="	",quote=FALSE,row.names=FALSE)
}
sessionInfo()
