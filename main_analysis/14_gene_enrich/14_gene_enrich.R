if (!exists(".shared_pipeline_call")) rm(list = ls())

gc()

getwd()

.libPaths()

if (!require("cli")) install.packages("cli")
library(Fast2TWAS)
library(data.table)

if (!exists("inputFile")) inputFile <- c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")

gtex_tissue_names <- fread("../reference/QTLEnrich/gtex_tissue_names.txt",header = F)$V1

for (trait in inputFile) {

  for (qtl_type in c("eQTL","sQTL")) {
    if (qtl_type=="eQTL") {
      QTLEnrich <- "../03_qtl_enrich_eqtl/result/GeneEnrich_Input_Files/GeneEnrich_input_best_eqtl_Tissue_"
    }else if(qtl_type=="sQTL"){
      QTLEnrich <- "../04_qtl_enrich_sqtl/result/GeneEnrich_Input_Files/GeneEnrich_input_best_sqtl_Tissue_"
    }

    for (tissue in gtex_tissue_names) {
      significant_genes <- paste0(QTLEnrich,tissue,"_",trait,"_",trait,"_significant.txt")
      null_genes <- paste0(QTLEnrich,tissue,"_",trait,"_",trait,"_null.txt")

      if (!file.exists(significant_genes) || !file.exists(null_genes) || length(significant_genes) != 1 || length(null_genes) != 1) {
        stop("Check QTLEnrich Info Info results")
      }
      gtf_file <- "../reference/GeneEnrich/gencode.v26.annotation.gtf"
      HLA_file <- "../reference/GeneEnrich/HLA_region_file_hg38.tsv"

      for (resource_name in c("GO_MF","GO_BP","GO_CC","KEGG","REACTOME","MGI","HALLMARK")) {
        gene_set_file <- paste0("../reference/GeneEnrich/",resource_name,"_genesets.txt")
        outname <- paste0(trait,"_",qtl_type,"_",tissue,"_",resource_name)

        outname_file <- list.files(path = "./result/", pattern = paste0("^GeneEnrich_", outname, ".*_log.txt$"), full.names = TRUE)
        if(length(outname_file) > 0 && file.exists(outname_file[1]))next

        significant_genes_header <- fread(significant_genes)
        null_genes_header <- fread(null_genes)

        setnames(significant_genes_header, old = c("eVariant", "eGene"), new = c("variant", "gene"))
        setnames(null_genes_header, old = c("eGene"), new = c("gene"))

        fwrite(significant_genes_header, "significant_genes_header.txt",sep="\t")
        fwrite(null_genes_header, "null_genes_header.txt",sep="\t")

        cli::cli_h1(paste0("Trait:",trait,"; Info qtl_type:",qtl_type,"; Tissue:",tissue,"; Info resource_name:",resource_name))

        GeneEnrich_run(
          significant_genes = "significant_genes_header.txt",
          null_genes = "null_genes_header.txt",
          gene_set_file = gene_set_file,
          min_permutations = 1000,
          max_permutations = 100000,
          gtf_file = gtf_file,
          resource_name = resource_name,
          is_restrict_genes_database = TRUE,
          is_HLA_remove = TRUE,
          HLA_file = HLA_file,
          gene_set_size_min = NULL,
          gene_set_size_max = NULL,
          is_fast_permute = TRUE,
          bed_remove = NULL,
          p_val_cutoff = NULL,
          q_val_cutoff = NULL,
          fdr_cutoff = NULL,
          genes_to_mark_table = NULL,
          mark_columns = NULL,
          gene_expression_levels = NULL,
          tissue_of_interest = NULL,
          is_development = FALSE,
          is_null_set = FALSE,
          is_sig_gs_cutoff = FALSE,
          is_help = FALSE,
          force = FALSE,
          outname = outname
        )

        unlink(c("significant_genes_header.txt","null_genes_header.txt"), recursive = T, force = T)

      }
    }
  }
}

sessionInfo()
