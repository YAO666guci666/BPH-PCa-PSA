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

write(x = "", file = "09_run_ecaviar.log", append = F)
for (trait in inputFile) {

  GenomicRiskLoci_dir <- list.files(
    path = paste0("../07_fuma_gwas/",trait,"/"),
    pattern = "^FUMA.*",
    full.names = TRUE,
    include.dirs = TRUE
  ) %>% .[file.info(.)$isdir]

  GenomicRiskLoci=fread(paste0(GenomicRiskLoci_dir,"/GenomicRiskLoci.txt"))

  GenomicRiskLoci$chr <- gsub("^23$|^x$", "X", as.character(GenomicRiskLoci$chr))

  for (snp in GenomicRiskLoci$rsID) {

    i=which(snp == GenomicRiskLoci$rsID)
    message("\n Trait:",trait,", Current",i,"GenomicRiskLoci,LeadSNPs:",GenomicRiskLoci[i,]$rsID)

    for (qtl_type in c("eQTL","sQTL")) {

      for (tissue in gtex_tissue_names) {

        message("\n Trait:",trait,", QTL type:",qtl_type,", Tissue:",tissue)

        trait_qtl_type_dir=paste0(trait,"/",qtl_type)
        trait_gwas_dir=paste0(trait,"/GWAS")
        trait_corr_dir=paste0(trait,"/corr")

        if (! dir.exists(paste0("result/",trait_qtl_type_dir))) {
          dir.create(paste0("result/",trait_qtl_type_dir), recursive = T)
        }

        qtl_data_file <- paste0("../08_prepare_ecaviar/",trait_qtl_type_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".qtl.txt")

        gwas_data_file <- paste0("../08_prepare_ecaviar/",trait_gwas_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".gwas.txt")

        ld_data_file <- paste0("../08_prepare_ecaviar/",trait_corr_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".corr.txt")

        if (!file.exists(qtl_data_file)  ||
            !file.exists(gwas_data_file) ||
            !file.exists(ld_data_file)   ||
            length(qtl_data_file) != 1   ||
            length(gwas_data_file) != 1  ||
            length(ld_data_file) != 1) {
          line=paste0("Trait:",trait,", Current lead SNP:",snp,", QTL type:",qtl_type,", Tissue:",tissue,", no overlapping SNPs")
          message(line)
          write(x = line, file = "09_run_ecaviar.log", append = T)
          next
        }

        caviar_eCAVIAR_run(gwas_z_file=gwas_data_file,
                           eqtl_z_file=qtl_data_file,
                           gwas_ld_file=ld_data_file,
                           eqtl_ld_file=ld_data_file,
                           rho_prob=0.95,
                           causal=2,
                           f=1,
                           outname=paste0(trait_qtl_type_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".eCAVIAR"))

      }
    }
  }
}

sessionInfo()
