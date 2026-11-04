if (!exists(".shared_pipeline_call")) rm(list = ls())

gc()

getwd()

.libPaths()

library(xQTLbiolinks)
library(Fast2TWAS)
library(data.table)

if (!exists("inputFile")) inputFile <- c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")

wind_size = 50000

ref_pop="g1000_eur"

gtex_tissue_names <- fread("../reference/QTLEnrich/gtex_tissue_names.txt",header = F)$V1

write(x = "", file = "13_process_fastenloc.log", append = F)
for (trait in inputFile) {

  if (!file.exists("./result")) {
    dir.create("./result")
  }

  stable6 <- fread(paste0("../07_fuma_gwas/stable6.txt"))

  stable6_filter <- subset(stable6, Trait == trait)

  gwas_data <- fread(paste0("../02_harmonize_gwas/",trait,"_hg38_qc.txt"))

  all_result <- data.frame()

  for (snp in unique(stable6_filter$rsID)) {

    for (qtl_type in c("eQTL","sQTL")) {

      for (tissue in gtex_tissue_names) {

        message("\n Trait:",trait,", Current lead SNP:",snp,", QTL type:",qtl_type,", Tissue:",tissue)

        if (file.exists(paste0("./result/stable_",snp,"_",trait,"_",qtl_type,"_",tissue,".fastenloc.txt")))next
        fastenloc_sig <- paste0("../12_run_fastenloc/result/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".fastenloc_result.enloc.sig.out")
        fastenloc_snp <- paste0("../12_run_fastenloc/result/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".fastenloc_result.enloc.snp.out")

        if (file.exists(fastenloc_sig) && file.exists(fastenloc_snp)) {

          fastenloc_sig_data <- read.table(fastenloc_sig,header = T)

          fastenloc_snp_data <- read.table(fastenloc_snp,header = T)  %>%
            tidyr::separate(
              col = SNP,
              into = c("gene_id_ENSG", "rsID_qtl"),
              sep = ":",
              remove = FALSE
            )

          query_with_retry <- function(rsID_qtl, max_retries = 99) {
            retries <- 0
            while (retries < max_retries) {
              tryCatch({
                varAsso_id <- xQTLquery_varId(rsID_qtl)
                return(varAsso_id)
              }, error = function(e) {
                message("Error encountered: ", e$message)
                retries <<- retries + 1
                message("Retrying... Attempt ", retries, " of ", max_retries)
                Sys.sleep(5)
              })
            }
            stop("Failed after ", max_retries, " attempts")
          }

          all_varAsso_id <- data.frame()
          for (rsID_qtl in fastenloc_snp_data$rsID_qtl) {

            varAsso_id <-try(query_with_retry(rsID_qtl), silent = TRUE)
            if (inherits(varAsso_id, "try-error")) {
              line=paste0("Error in querying variant ID for SNP: ", rsID_qtl)
              write(x = line, file = "13_process_fastenloc.log", append = T)
              message("Error in querying variant ID for SNP: ", rsID_qtl)
              next
            }
            all_varAsso_id <- rbind(all_varAsso_id, varAsso_id)
          }

          if (nrow(all_varAsso_id)==0) next
          all_varAsso_id <- subset(all_varAsso_id, select=c(variantId, snpId,ref,alt))

          if (qtl_type=="eQTL") {
            qtl_data_signif_variant_file <- paste0("../reference/QTLEnrich/GTEx_Analysis_v8_",qtl_type,"/",tissue,".v8.signif_variant_gene_pairs.txt.gz")
            qtl_data_signif_variant <- fread(qtl_data_signif_variant_file)

            qtl_data_signif_variant_fileter <- subset(qtl_data_signif_variant,variant_id %in% all_varAsso_id$variantId,
                                                      select=c(variant_id,gene_id,pval_nominal,slope,slope_se))
            qtl_data_signif_variant_fileter$Gene_id_intron_cluster <- ""
          }else if(qtl_type=="sQTL"){
            qtl_data_signif_variant_file <- paste0("../reference/QTLEnrich/GTEx_Analysis_v8_",qtl_type,"/",tissue,".v8.sqtl_signifpairs.txt.gz")
            qtl_data_signif_variant <- fread(qtl_data_signif_variant_file)
            setnames(qtl_data_signif_variant,"phenotype_id","Gene_id_intron_cluster")
            qtl_data_signif_variant <- qtl_data_signif_variant %>% tidyr::separate(col = Gene_id_intron_cluster,
                                                                                   into = c("chr", "start", "end", "clu","gene_id"),
                                                                                   sep = ":",
                                                                                   remove = FALSE
            )

            qtl_data_signif_variant_fileter <- subset(qtl_data_signif_variant,variant_id %in% all_varAsso_id$variantId,
                                                      select=c(variant_id,Gene_id_intron_cluster,gene_id,pval_nominal,slope,slope_se))

          }

          all_varAsso_id_qtl_data_signif_variant_fileter <- merge(qtl_data_signif_variant_fileter, all_varAsso_id, by.x="variant_id", by.y="variantId")

          fastenloc_result <- merge(fastenloc_snp_data, all_varAsso_id_qtl_data_signif_variant_fileter, by.x="rsID_qtl", by.y="snpId")

          fastenloc_result$Tissue <- tissue

          GTEXv8_ensTOsymble <- GTEXv8_ensTOsymble[,c("gene_id","symble")]

          fastenloc_result$gene_id <- gsub("\\..*", "", fastenloc_result$gene_id)
          GTEXv8_ensTOsymble_fastenloc_result <- merge(fastenloc_result, GTEXv8_ensTOsymble, by.x="gene_id", by.y="gene_id")

          gwas_data_filter <- subset(gwas_data, SNP %in% unique(fastenloc_result$rsID_qtl),select = c(SNP,effect_allele,other_allele,pval,beta,se))

          GTEXv8_ensTOsymble_fastenloc_result_final <- merge(GTEXv8_ensTOsymble_fastenloc_result, gwas_data_filter, by.x="rsID_qtl", by.y="SNP")

          for (i in 1:nrow(GTEXv8_ensTOsymble_fastenloc_result_final)) {

            if (GTEXv8_ensTOsymble_fastenloc_result_final$effect_allele[i]==GTEXv8_ensTOsymble_fastenloc_result_final$alt[i] &
                GTEXv8_ensTOsymble_fastenloc_result_final$other_allele[i]==GTEXv8_ensTOsymble_fastenloc_result_final$ref[i]) {

              if (GTEXv8_ensTOsymble_fastenloc_result_final$slope[i]>0 & GTEXv8_ensTOsymble_fastenloc_result_final$beta[i]>0) {
                GTEXv8_ensTOsymble_fastenloc_result_final$Qtl_gwas_relative_effect[i] <- "SAME"
              } else if (GTEXv8_ensTOsymble_fastenloc_result_final$slope[i]<0 & GTEXv8_ensTOsymble_fastenloc_result_final$beta[i]<0) {
                GTEXv8_ensTOsymble_fastenloc_result_final$Qtl_gwas_relative_effect[i] <- "SAME"
              } else {
                GTEXv8_ensTOsymble_fastenloc_result_final$Qtl_gwas_relative_effect[i] <- "OPPOSITE"
              }
            }else if(GTEXv8_ensTOsymble_fastenloc_result_final$effect_allele[i]==GTEXv8_ensTOsymble_fastenloc_result_final$ref[i]  &
                     GTEXv8_ensTOsymble_fastenloc_result_final$other_allele[i]==GTEXv8_ensTOsymble_fastenloc_result_final$alt[i]) {

              if (GTEXv8_ensTOsymble_fastenloc_result_final$slope[i]>0 & GTEXv8_ensTOsymble_fastenloc_result_final$beta[i]>0) {
                GTEXv8_ensTOsymble_fastenloc_result_final$Qtl_gwas_relative_effect[i] <- "OPPOSITE"
              } else if (GTEXv8_ensTOsymble_fastenloc_result_final$slope[i]<0 & GTEXv8_ensTOsymble_fastenloc_result_final$beta[i]<0) {
                GTEXv8_ensTOsymble_fastenloc_result_final$Qtl_gwas_relative_effect[i] <- "OPPOSITE"
              } else {
                GTEXv8_ensTOsymble_fastenloc_result_final$Qtl_gwas_relative_effect[i] <- "SAME"
              }
            }else{
              GTEXv8_ensTOsymble_fastenloc_result_final$Qtl_gwas_relative_effect[i] <- "NA"
            }
          }

          stable6_filter_snp <- subset(stable6_filter, rsID==snp, select=-c(hg19,CHR,N))
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "variant_id", "qtl_Variant_id")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "symble", "gene_symble")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "pval_nominal", "QTL_pval_nominal")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "pval", "QTL_gwas_p_value")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "alt", "qtl_variant_ALT_(effect_allele)")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "ref", "qtl_variant_REF_(other_allele)")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "slope", "QTL_slope_relative_ALT")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "slope_se", "QTL_se")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "beta", "GWAS_effect_relative_ALT")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "se", "GWAS_se_relative_ALT")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "effect_allele", "GWAS_variant_ALT_(effect_allele)")
          setnames(GTEXv8_ensTOsymble_fastenloc_result_final, "other_allele", "GWAS_variant_REF_(other_allele)")

          GTEXv8_ensTOsymble_fastenloc_result_final <- merge(GTEXv8_ensTOsymble_fastenloc_result_final,fastenloc_sig_data)

          GTEXv8_ensTOsymble_fastenloc_result_final <- subset(GTEXv8_ensTOsymble_fastenloc_result_final,select = -c(Signal,SNP,gene_id_ENSG,PIP_qtl,PIP_gwas_marginal,SCP,CPIP_qtl,CPIP_gwas_marginal,CPIP_gwas_qtl_prior,LCP))

          setnames(stable6_filter_snp, "rsID", "rsID_GWAS")
          setnames(stable6_filter_snp, "GWAS_locus", "GWAS_variant")
          setnames(stable6_filter_snp, "beta", "GWAS_effect")
          setnames(stable6_filter_snp, "se", "GWAS_se")
          setnames(stable6_filter_snp, "pval", "GWAS_p_value")
          setnames(stable6_filter_snp, "effect_allele", "GWAS_effect_allele")
          setnames(stable6_filter_snp, "other_allele", "GWAS_other_allele")

          stable6_filter_snp$startbp <- stable6_filter_snp$hg38 - wind_size
          stable6_filter_snp$endbp <- stable6_filter_snp$hg38 + wind_size

          stable6_filter_snp_repeated  <- stable6_filter_snp[rep(1,nrow(GTEXv8_ensTOsymble_fastenloc_result_final)),]
          GTEXv8_ensTOsymble_fastenloc_result_final_stable6 <- cbind(stable6_filter_snp_repeated,GTEXv8_ensTOsymble_fastenloc_result_final)

          GTEXv8_ensTOsymble_fastenloc_result_final_stable6$QTL_Type <- qtl_type
          all_result <- rbind(all_result,GTEXv8_ensTOsymble_fastenloc_result_final_stable6)

          fwrite(GTEXv8_ensTOsymble_fastenloc_result_final_stable6,
                 file=paste0("./result/stable_",snp,"_",trait,"_",qtl_type,"_",tissue,".fastenloc.txt"),
                 sep = "\t")

          GTEXv8_ensTOsymble_fastenloc_result_final_stable6 <- data.frame()

        } else {
          line=paste0("File not found:", fastenloc_sig, "or", fastenloc_snp)
          write(x = line, file = "13_process_fastenloc.log", append = T)
        }
      }
    }
  }

  fwrite(all_result,
         file=paste0("stable_",trait,".fastenloc.txt"),
         sep = "\t")

}

for (trait in inputFile) {

  data_file <- list.files(path = "./result/", pattern = paste0("^stable_.*.", trait, ".*.fastenloc.txt$"), full.names = TRUE)

  data <- lapply(data_file, fread)

  data <- do.call(rbind,data)

  fwrite(data,
         file=paste0("stable_",trait,".fastenloc.txt"),
         sep = "\t")

  data_filter <- subset(data,RCP>0.1)

  fwrite(data_filter,
         file=paste0("stable_",trait,".fastenloc.sig.txt"),
         sep = "\t")

}

sessionInfo()
