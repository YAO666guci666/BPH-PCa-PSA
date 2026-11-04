if (!exists(".shared_pipeline_call")) rm(list = ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)
library(data.table)
library(dplyr)
library(ggplot2)

if (!exists("inputFile")) inputFile <- c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")

gtex_tissue_names <- fread("../reference/QTLEnrich/gtex_tissue_names.txt",header = F)$V1

width = 10
height = 8

for (trait in inputFile) {

  all_GeneEnrich_result <- data.frame()

  for (qtl_type in c("eQTL","sQTL")) {

    for (tissue in gtex_tissue_names) {

      for (resource_name in c("GO_MF","GO_BP","GO_CC","KEGG","REACTOME","MGI","HALLMARK")) {
        outname <- paste0(trait,"_",qtl_type,"_",tissue,"_",resource_name)

        outname_file <- list.files(path = "../14_gene_enrich/result/", pattern = paste0("^GeneEnrich_gs_results_", outname, ".*.tsv$"), full.names = TRUE)
        if(length(outname_file) > 0 && file.exists(outname_file[1])){

          cli::cli_h2(paste0(outname_file,"already exists, reading it directly"))
          GeneEnrich_result <- fread(outname_file[1])

          GeneEnrich_result$Trait <- trait

          GeneEnrich_result$Tissue <- tissue

          GeneEnrich_result$QTL_type <- qtl_type

          all_GeneEnrich_result <- rbind(all_GeneEnrich_result,GeneEnrich_result)

        }
      }
    }
  }

  fwrite(all_GeneEnrich_result,
         file=paste0("stable_",trait,".GSEA_QTLEnrich.txt"),
         sep = "\t")

}

sessionInfo()
