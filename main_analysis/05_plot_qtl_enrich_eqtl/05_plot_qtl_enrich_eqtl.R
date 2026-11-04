rm(list=ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)
library(data.table)

inputFile=c("finngen_BPH","GCST90476277_PSA","ieu-b-4809_PC")

tissue_data <-fread("../reference/QTLEnrich/gtex_tissue_names.txt",header = F)
for (trait in inputFile) {
  for (tissue in tissue_data$V1) {

  message("\n Trait:",trait,"; Current",which(tissue == tissue_data$V1),"tissue:",tissue)

  key <- paste0(trait,"_",tissue)

  qtl_file <- list.files(path = "../03_qtl_enrich_eqtl/result/output_files/",
                         pattern = paste0("^QTLEnrich_significant.*",key),
                         full.names = T)
  pvalue_null_matrix  <- list.files(path = "../03_qtl_enrich_eqtl/result/output_files/",
                                    pattern = paste0("^QTLEnrich_matched_null.*",key),
                                    full.names = T)

  if (!file.exists(qtl_file) || !file.exists(pvalue_null_matrix) || length(qtl_file) != 1 || length(pvalue_null_matrix) != 1) {
    stop("Check QTLEnrich Info Info results")
  }

  QTLEnrich_qqplot_run(
    qtl_file=qtl_file,
    gwas_file=paste0("../02_harmonize_gwas/",trait,"_hg38_qc_QTLEnrich.txt"),
    pvalue_null_matrix=pvalue_null_matrix,
    is_gwas_background=T,
    is_matched_null_background=T,
    num_permutations=100,
    is_confidence_interval=T,
    upper_bound_ci=97.5,
    lower_bound_ci=2.5,
    upper_bound_clip=0,
    trait_label=trait,
    qtl_type="eQTLs",
    tissue_label=tissue,
    is_show_title=TRUE,
    qqplot_title=paste0(trait," GWAS"),
    flag_options=NULL,
    x_inches=12,
    y_inches=10,
    dpi=1000,
    exp_label=trait,
    is_help=F,
    force=F)

  }
}

sessionInfo()
