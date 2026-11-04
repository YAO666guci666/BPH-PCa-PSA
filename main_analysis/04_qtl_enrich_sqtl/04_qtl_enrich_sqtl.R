if (!exists(".shared_pipeline_call")) rm(list = ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)

if (!exists("inputFile")) inputFile <- c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")

for (trait in inputFile) {
  for (qtl_type in c("best_sqtl")) {
  QTLEnrich_run(gwas_file=paste0("../02_harmonize_gwas/",trait,"_hg38_qc_QTLEnrich.txt"),
                qtl_directory="../reference/QTLEnrich/GTEx_Analysis_v8_sQTL/",
                file_name=".v8.sgenes.txt.gz",
                qtl_type=qtl_type,
                confounders_table="../reference/QTLEnrich/GTEx_v8_Confounders_Table_49Tiss_Global.txt.gz",
                null_table="../reference/QTLEnrich/GTEx_v8_Null_Table_49Tiss_Global.txt.gz",
                is_subset_genes=TRUE,
                is_compute_tss_distance=TRUE,
                is_keep_gene_id_suffix=TRUE,
                is_keep_null_variants=TRUE,
                is_keep_pvalue_matrix=TRUE,
                gencode_file="../reference/QTLEnrich/GENCODE_26_df.txt",
                null_option=NULL,
                genome_build=NULL,
                trait_name=trait,
                qtl_q_value=NULL,
                gwas_p_value=NULL,
                independent_ranking=NULL,
                lambda_factor=NULL,
                lower_bound_permutations=NULL,
                upper_bound_permutations=NULL,
                num_quantiles=NULL,
                tissue_names=NULL,
                is_subset_tissues=FALSE,
                exp_label=trait,
                is_GeneEnrich_input=TRUE,
                eGenes_directory=NULL,
                eGenes_file_name=NULL,
                flag_options=NULL,
                is_help=FALSE,
                force=FALSE)
  }
}

sessionInfo()
