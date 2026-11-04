rm(list=ls())
gc()
library(Fast2TWAS)
library(data.table)

inputFile="GCST90476319.tsv.gz"
trait="GCST90476319_PSA"
data=fread(paste0("../rawdata/",inputFile))
print(names(data))
data=gwas2format(gwas_dat=data,
  chr_col="chromosome",pos_col="base_pair_location",snp_col="rsid",
  effect_allele_col="effect_allele",other_allele_col="other_allele",
  eaf_col="effect_allele_frequency",beta_col="beta",se_col="standard_error",
  pval_col="p_value",log_pval=FALSE,min_pval=1e-200,samplesize=254502)

data[, c("effect_allele","other_allele") := .(other_allele,effect_allele)]
data[, eaf := 1 - eaf]
data[, beta := -beta]
fwrite(data,file=paste0(trait,".txt"),sep="	")
fwrite(data,file="PSA_continuous.txt",sep="\t")
sessionInfo()
