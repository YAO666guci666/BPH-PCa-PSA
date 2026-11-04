rm(list=ls())
gc()
library(Fast2TWAS)
library(data.table)

inputFile="finngen_BPH.gz"
trait="finngen_BPH"
N=194710
data=finngenClean(finngenfile = paste0("../rawdata/",inputFile),
                  samplesize = N,
                  savefile = paste0(trait,".txt"))
data$ncase <- 41137
data$ncontrol <- 153573
fwrite(data,file=paste0(trait,".txt"),sep="	")
fwrite(data,file="BPH.txt",sep="\t")

inputFile="ieu-b-4809.vcf.gz"
trait="ieu-b-4809_PC"
data=read_ieu_data(ieu_id=paste0("../rawdata/",inputFile),
                   download=FALSE,maf=0.01,pvalue=NULL,
                   is_fast2twas_format=TRUE)
data$ncase <- 9132
data$ncontrol <- 173493
fwrite(data,file=paste0(trait,".txt"),sep="	")
fwrite(data,file="PCA.txt",sep="\t")

inputFile="GCST90476277_PSA.txt.gz"
trait="GCST90476277_PSA"
data=fread(paste0("../rawdata/",inputFile))
print(names(data))
data=gwas2format(gwas_dat=data,
  chr_col="chromosome",pos_col="base_pair_location",snp_col="rsid",
  effect_allele_col="effect_allele",other_allele_col="other_allele",
  eaf_col="effect_allele_frequency",beta_col="beta",se_col="standard_error",
  pval_col="p_value",log_pval=FALSE,min_pval=1e-200,samplesize=401453)
data$ncase <- 54381
data$ncontrol <- 347072
fwrite(data,file=paste0(trait,".txt"),sep="	")
fwrite(data,file="PSA.txt",sep="\t")
sessionInfo()
