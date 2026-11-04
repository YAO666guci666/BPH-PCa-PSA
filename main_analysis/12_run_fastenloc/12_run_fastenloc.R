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

thread = 8

write(x = "", file = "12_run_fastenloc.log", append = F)
for (trait in inputFile) {

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

        trait_qtl_type_dir=paste0("../11_prepare_fastenloc/",trait,"/")

        gwas_qtl_all_file <- paste0(trait_qtl_type_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".fastenloc.txt")
        if (!file.exists(gwas_qtl_all_file)) {
          line=paste0("Trait:",trait,", Current lead SNP:",snp,", QTL type:",qtl_type,", Tissue:",tissue,", does not exist file")
          message(line)
          write(x = line, file = "12_run_fastenloc.log", append = T)
          next
        }
        gwas_qtl_all_data <- fread(gwas_qtl_all_file)
        total_variant <- nrow(gwas_qtl_all_data)
        outname_file <- paste0(snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".fastenloc_result")

        fastenloc_coloc_run(summary = gwas_qtl_all_file,
                            total_variant = total_variant,
                            thread = thread,
                            is_coloc_default_prior = T,
                            outname = outname_file)

      }
    }
  }
}

sessionInfo()
