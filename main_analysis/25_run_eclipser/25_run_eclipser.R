rm(list=ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)
library(data.table)

sc_path <- "../23_prepare_eclipser/PC_train.csv"

tissues <- unique(fread(sc_path)$"tissue")

clumping_files <- list.files(path = "../24_clump_eclipser/result/",
                             pattern = "*consolidated_clumping.tsv",
                             recursive = TRUE,
                             full.names = TRUE)

for (clumping_file in clumping_files) {
  for (tissue in tissues) {
    if (basename(clumping_file)=="Null_consolidated_clumping.tsv")next
    trait <- gsub("_consolidated_clumping.tsv", "", basename(clumping_file))

    ECLIPSER_run(sc_path=sc_path,
                 clumps_file=clumping_file,
                 background_file="../24_clump_eclipser/result/Null_consolidated_clumping.tsv" ,
                 trait=trait,
                 tissue=tissue,
                 num_samples=100,000,
                 method="Bayesian",
                 cutoff1=0.50,
                 cutoff2=0.1,
                 cutoff3=0.05)
  }
}

sessionInfo()
