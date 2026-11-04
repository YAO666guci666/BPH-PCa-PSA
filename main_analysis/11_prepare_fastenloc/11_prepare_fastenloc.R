if (!exists(".shared_pipeline_call")) rm(list = ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)
library(data.table)

if (!exists("inputFile")) inputFile <- c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")

wind_size = 50000

ref_pop="g1000_eur"

gtex_tissue_names <- fread("../reference/QTLEnrich/gtex_tissue_names.txt",header = F)$V1

write(x = "", file = "11_prepare_fastenloc.log", append = F)
for (trait in inputFile) {

  trait_qtl_type_dir=paste0(trait)
  if (! dir.exists(trait_qtl_type_dir)) {
    dir.create(trait_qtl_type_dir, recursive = T)
  }

  GenomicRiskLoci_dir <- list.files(
    path = paste0("../07_fuma_gwas/",trait,"/"),
    pattern = "^FUMA.*",
    full.names = TRUE,
    include.dirs = TRUE
  ) %>% .[file.info(.)$isdir]

  GenomicRiskLoci=fread(paste0(GenomicRiskLoci_dir,"/GenomicRiskLoci.txt"))

  for (snp in GenomicRiskLoci$rsID) {

    i=which(snp == GenomicRiskLoci$rsID)
    message("\n Trait:",trait,", Current",i,"GenomicRiskLoci,LeadSNPs:",GenomicRiskLoci[i,]$rsID)

    for (qtl_type in c("eQTL","sQTL")) {

      for (tissue in gtex_tissue_names) {

        message("\n Trait:",trait,", QTL type:",qtl_type,", Tissue:",tissue,", Current lead SNP:",snp)

        trait_gwas_dir=paste0("../08_prepare_ecaviar/",trait,"/GWAS")

        gwas_qtl_all_file <- paste0(trait_gwas_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".all.txt")
        if (!file.exists(gwas_qtl_all_file)) {
          line=paste0("Trait:",trait,", Current lead SNP:",snp,", QTL type:",qtl_type,", Tissue:",tissue,", does not exist file")
          message(line)
          write(x = line, file = "11_prepare_fastenloc.log", append = T)
          next
        }

        markers_QTL_gwas_final_SNP <- fread(gwas_qtl_all_file)

        for (i in 1:nrow(markers_QTL_gwas_final_SNP)) {

            if(markers_QTL_gwas_final_SNP$effect_allele[i]==markers_QTL_gwas_final_SNP$ref[i]  &
                   markers_QTL_gwas_final_SNP$other_allele[i]==markers_QTL_gwas_final_SNP$alt[i]) {
              markers_QTL_gwas_final_SNP$slope[i] <- -markers_QTL_gwas_final_SNP$slope[i]
            }
        }

        if (qtl_type == "sQTL") {

          markers_QTL_gwas_final_SNP <- markers_QTL_gwas_final_SNP %>%
            tidyr::separate(
              col = phenotype_id,
              into = c("chr", "start", "end", "cul","gene_id"),
              sep = ":",
              remove = FALSE
            )
        }

        markers_QTL_gwas_final_SNP <- markers_QTL_gwas_final_SNP[,c("gene_id","SNP","slope","slope_se","beta","se")]

        fwrite(markers_QTL_gwas_final_SNP,
               file=paste0(trait_qtl_type_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".fastenloc.txt"),
               sep = "\t",
               col.names = F)
      }
    }
  }
}

sessionInfo()
