if (!exists(".shared_pipeline_call")) rm(list = ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)
library(data.table)

if (!exists("inputFile")) inputFile <- c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")

gtex_tissue_names <- fread("../reference/QTLEnrich/gtex_tissue_names.txt",header = F)$V1

for (trait in inputFile) {

  all_base_data <- data.frame()

  stable6 <- fread(paste0("../07_fuma_gwas/stable6.txt"))

  stable6_filter <- subset(stable6, Trait == trait)

  all_eCAVIAR_result <- fread(paste0("../10_process_ecaviar/stable_eCAVIAR_",trait,".txt"))
  all_fastENLOC_result <- fread(paste0("../13_process_fastenloc/stable_",trait,".fastenloc.txt"))

  FUMA_path <- list.files(
    path = paste0("../07_fuma_gwas/",trait,"/"),
    pattern = "^FUMA.*",
    full.names = TRUE,
    include.dirs = TRUE
  ) %>% .[file.info(.)$isdir]

  dat_snps=data.table::fread(paste0(FUMA_path,"/snps.txt"))[,c("rsID","nearestGene")]

  for (snp in unique(stable6_filter$rsID)) {

    eCAVIAR_result_filter_snp <- subset(all_eCAVIAR_result, rsID_GWAS == snp)

    eCAVIAR_result_filter_snp_eQTL <- subset(eCAVIAR_result_filter_snp, QTL_Type == "eQTL")
    eQTL_tissue_geng_pairs_tested_eCAVIAR <- nrow(eCAVIAR_result_filter_snp_eQTL)
    eGenes_tested_eCAVIAR <- length(unique(eCAVIAR_result_filter_snp_eQTL$gene_symble))

    significant_eCAVIAR_result_filter_snp_eQTL <- subset(eCAVIAR_result_filter_snp_eQTL, CLPP >0.01)
    significant_eQTL_tissue_geng_pairs_tested_eCAVIAR <- nrow(significant_eCAVIAR_result_filter_snp_eQTL)
    significant_eGenes_eCAVIAR <- length(unique(significant_eCAVIAR_result_filter_snp_eQTL$gene_symble))
    List_of_colocalizing_eGenes_eCAVIAR <- paste(unique(significant_eCAVIAR_result_filter_snp_eQTL$gene_symble),collapse = ", ")
    if (List_of_colocalizing_eGenes_eCAVIAR=="") {
      List_of_colocalizing_eGenes_eCAVIAR=NA
    }

    eCAVIAR_result_filter_snp_sQTL <- subset(eCAVIAR_result_filter_snp, QTL_Type == "sQTL")
    sQTL_tissue_geng_pairs_tested_eCAVIAR <- nrow(eCAVIAR_result_filter_snp_sQTL)
    sGenes_tested_eCAVIAR <- length(unique(eCAVIAR_result_filter_snp_sQTL$gene_symble))

    significant_eCAVIAR_result_filter_snp_sQTL <- subset(eCAVIAR_result_filter_snp_sQTL, CLPP >0.01)
    significant_sQTL_tissue_geng_pairs_tested_eCAVIAR <- nrow(significant_eCAVIAR_result_filter_snp_sQTL)
    significant_sGenes_eCAVIAR <- length(unique(significant_eCAVIAR_result_filter_snp_sQTL$gene_symble))
    List_of_colocalizing_sGenes_eCAVIAR <- paste(unique(significant_eCAVIAR_result_filter_snp_sQTL$gene_symble),collapse = ", ")
    if (List_of_colocalizing_sGenes_eCAVIAR=="") {
      List_of_colocalizing_sGenes_eCAVIAR=NA
    }
    significant_sGenes_by_intron_cluster_eCAVIAR <- length(unique(significant_eCAVIAR_result_filter_snp_sQTL$Gene_id_intron_cluster))
    List_of_colocalizing_sGenes_by_intron_cluster_eCAVIAR <- paste(unique(significant_eCAVIAR_result_filter_snp_sQTL$Gene_id_intron_cluster),collapse = ", ")
    if (List_of_colocalizing_sGenes_by_intron_cluster_eCAVIAR=="") {
      List_of_colocalizing_sGenes_by_intron_cluster_eCAVIAR=NA
    }

    eGenes_sGenes_tested_eCAVIAR <- length(unique(eCAVIAR_result_filter_snp$gene_symble))
    significantly_colocalizing_eGenes_sGenes_eCAVIAR <- length(unique(c(significant_eCAVIAR_result_filter_snp_eQTL$gene_symble,
                                                                        significant_eCAVIAR_result_filter_snp_sQTL$gene_symble)
                                                                      ))
    if (is.na(List_of_colocalizing_sGenes_eCAVIAR)) {
      sGenes_eCAVIAR <- NULL
    } else {
      sGenes_eCAVIAR <- paste0(unique(significant_eCAVIAR_result_filter_snp_sQTL$gene_symble), "(s)",collapse = ", ")
    }

    if (is.na(List_of_colocalizing_eGenes_eCAVIAR)) {
      eGenes_eCAVIAR <- NULL
    } else {
      eGenes_eCAVIAR <- paste0(unique(significant_eCAVIAR_result_filter_snp_eQTL$gene_symble), "(e)",collapse = ", ")
    }

    List_of_significantly_colocalizing_eGenes_sGenes_eCAVIAR <- paste(c(eGenes_eCAVIAR,sGenes_eCAVIAR),collapse = ", ")

    if(length(significant_eCAVIAR_result_filter_snp_eQTL$Tissue)==0){
      tissues_eGenes_eCAVIAR <- NULL
    }else{
      tissues_eGenes_eCAVIAR <- paste0(unique(significant_eCAVIAR_result_filter_snp_eQTL$Tissue), "(e)",collapse = ", ")
    }
    if(length(significant_eCAVIAR_result_filter_snp_sQTL$Tissue)==0){
      tissues_sGenes_eCAVIAR <- NULL
    }else{
      tissues_sGenes_eCAVIAR <- paste0(unique(significant_eCAVIAR_result_filter_snp_sQTL$Tissue), "(s)",collapse = ", ")
    }

    List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_eCAVIAR <- paste(c(tissues_eGenes_eCAVIAR,tissues_sGenes_eCAVIAR),collapse = ", ")

    significantly_fastENLOC_result_filter_snp <- subset(all_fastENLOC_result, rsID_GWAS == snp & RCP > 0.1)

    significantly_fastENLOC_result_filter_snp_eQTL <- subset(significantly_fastENLOC_result_filter_snp, QTL_Type == "eQTL")
    significantly_eQTL_tissue_geng_pairs_tested_fastENLOC <- nrow(significantly_fastENLOC_result_filter_snp_eQTL)
    significantly_eGenes_tested_fastENLOC <- length(unique(significantly_fastENLOC_result_filter_snp_eQTL$gene_symble))
    List_of_significantly_colocalizing_eGenes_fastENLOC <- paste(unique(significantly_fastENLOC_result_filter_snp_eQTL$gene_symble),collapse = ", ")
    if (List_of_significantly_colocalizing_eGenes_fastENLOC=="") {
      List_of_significantly_colocalizing_eGenes_fastENLOC=NA
    }

    significantly_fastENLOC_result_filter_snp_sQTL <- subset(significantly_fastENLOC_result_filter_snp, QTL_Type == "sQTL")
    significantly_sQTL_tissue_geng_pairs_tested_fastENLOC <- nrow(significantly_fastENLOC_result_filter_snp_sQTL)
    significantly_sGenes_tested_fastENLOC <- length(unique(significantly_fastENLOC_result_filter_snp_sQTL$gene_symble))
    List_of_significantly_colocalizing_sGenes_fastENLOC <- paste(unique(significantly_fastENLOC_result_filter_snp_sQTL$gene_symble),collapse = ", ")
    if (List_of_significantly_colocalizing_sGenes_fastENLOC=="") {
      List_of_significantly_colocalizing_sGenes_fastENLOC=NA
    }

    significant_sGenes_by_intron_cluster_fastENLOC <- length(unique(significantly_fastENLOC_result_filter_snp_sQTL$Gene_id_intron_cluster))
    List_of_significantly_colocalizing_sGenes_by_intron_cluster_fastENLOC <- paste(unique(significantly_fastENLOC_result_filter_snp_sQTL$Gene_id_intron_cluster),collapse = ", ")
    if (List_of_significantly_colocalizing_sGenes_by_intron_cluster_fastENLOC=="") {
      List_of_significantly_colocalizing_sGenes_by_intron_cluster_fastENLOC=NA
    }

    significantly_colocalizing_eGenes_sGenes_fastENLOC <- length(unique(significantly_fastENLOC_result_filter_snp$gene_symble))
    if (is.na(List_of_significantly_colocalizing_sGenes_fastENLOC)) {
      significantly_sGenes_fastENLOC <- NULL
    } else {
      significantly_sGenes_fastENLOC <- paste0(unique(significantly_fastENLOC_result_filter_snp_sQTL$gene_symble), "(s)",collapse = ", ")
    }

    if (is.na(List_of_significantly_colocalizing_eGenes_fastENLOC)) {
      significantly_eGenes_fastENLOC <- NULL
    } else {
      significantly_eGenes_fastENLOC <- paste0(unique(significantly_fastENLOC_result_filter_snp_eQTL$gene_symble), "(e)",collapse = ", ")
    }

    List_of_significantly_colocalizing_eGenes_sGenes_fastENLOC <- paste(c(significantly_eGenes_fastENLOC,significantly_sGenes_fastENLOC),collapse = ", ")

    if(length(significantly_fastENLOC_result_filter_snp_eQTL$Tissue)==0){
      tissues_eGenes_fastENLOC <- NULL
    }else{
      tissues_eGenes_fastENLOC <- paste0(unique(significantly_fastENLOC_result_filter_snp_eQTL$Tissue), "(e)",collapse = ", ")
    }
    if(length(significantly_fastENLOC_result_filter_snp_sQTL$Tissue)==0){
      tissues_sGenes_fastENLOC <- NULL
    }else{
      tissues_sGenes_fastENLOC <- paste0(unique(significantly_fastENLOC_result_filter_snp_sQTL$Tissue), "(s)",collapse = ", ")
    }

    List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_fastENLOC <- paste(c(tissues_eGenes_fastENLOC,tissues_sGenes_fastENLOC),collapse = ", ")

    significant_eQTL_tissue_geng_pairs_tested_eCAVIAR_union <- significant_eCAVIAR_result_filter_snp_eQTL[,c("rsID_qtl","Tissue","gene_symble")]
    significant_eQTL_tissue_geng_pairs_tested_fastENLOC_union <- significantly_fastENLOC_result_filter_snp_eQTL[,c("rsID_qtl","Tissue","gene_symble")]

    union_eQTL_eCAVIAR_fastENLOC_union <- rbind(significant_eQTL_tissue_geng_pairs_tested_eCAVIAR_union,significant_eQTL_tissue_geng_pairs_tested_fastENLOC_union)
    significant_eQTL_tissue_geng_pairs_union <- unique(union_eQTL_eCAVIAR_fastENLOC_union)
    significant_eGenes_tested_union <- length(unique(significant_eQTL_tissue_geng_pairs_union$gene_symble))
    List_of_significant_colocalizing_eGenes_union <- paste(unique(c(significant_eQTL_tissue_geng_pairs_union$gene_symble)),collapse = ", ")
    if (List_of_significant_colocalizing_eGenes_union=="") {
      List_of_significant_colocalizing_eGenes_union=NULL
    }

    significant_sQTL_tissue_geng_pairs_tested_eCAVIAR_union <- significant_eCAVIAR_result_filter_snp_sQTL[,c("rsID_qtl","Tissue","gene_symble")]
    significant_sQTL_tissue_geng_pairs_tested_fastENLOC_union <- significantly_fastENLOC_result_filter_snp_sQTL[,c("rsID_qtl","Tissue","gene_symble")]

    union_sQTL_eCAVIAR_fastENLOC_union <- rbind(significant_sQTL_tissue_geng_pairs_tested_eCAVIAR_union,significant_sQTL_tissue_geng_pairs_tested_fastENLOC_union)
    significant_sQTL_tissue_geng_pairs_union <- unique(union_sQTL_eCAVIAR_fastENLOC_union)
    significant_sGenes_tested_union <- length(unique(significant_sQTL_tissue_geng_pairs_union$gene_symble))
    List_of_significant_colocalizing_sGenes_union <- paste(unique(significant_sQTL_tissue_geng_pairs_union$gene_symble),collapse = ", ")
    if (List_of_significant_colocalizing_sGenes_union=="") {
      List_of_significant_colocalizing_sGenes_union=NULL
    }

    List_of_significant_colocalizing_sGenes_by_intron_cluster_union <- paste(c(unique(significantly_fastENLOC_result_filter_snp_sQTL$Gene_id_intron_cluster,
                                                                                      significant_eCAVIAR_result_filter_snp_sQTL$Gene_id_intron_cluster)),
                                                                             collapse = ", ")
    if (List_of_significant_colocalizing_sGenes_by_intron_cluster_union=="") {
      List_of_significant_colocalizing_sGenes_by_intron_cluster_union=NA
    }

    significant_sGenes_by_intron_cluster_union <- length(unique(significantly_fastENLOC_result_filter_snp_sQTL$Gene_id_intron_cluster,
                                                                significant_eCAVIAR_result_filter_snp_sQTL$Gene_id_intron_cluster))

    significantly_colocalizing_eGenes_sGenes_union <- length(unique(c(significant_eCAVIAR_result_filter_snp_eQTL$gene_symble,
                                                                      significant_eCAVIAR_result_filter_snp_sQTL$gene_symble,
                                                                      significantly_fastENLOC_result_filter_snp_eQTL$gene_symble,
                                                                      significantly_fastENLOC_result_filter_snp_sQTL$gene_symble)
                                                                    ))
    List_of_significantly_colocalizing_eGenes_sGenes_union <- paste(unique(c(paste0(unlist(strsplit(as.character(List_of_significant_colocalizing_eGenes_union), ", "))),
                                                                             paste0(unlist(strsplit(as.character(List_of_significant_colocalizing_sGenes_union), ", "))))),
                                                                    collapse = ", ")

    List_of_significantly_colocalizing_tissues_from_eGenes_sGenes_union <- paste(unique(c(unlist(strsplit(as.character(List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_fastENLOC), ", ")),
                                                                                          unlist(strsplit(as.character(List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_eCAVIAR), ", ")))),
                                                                                 collapse = ", ")

    significant_eGenes_intersection <- length(intersect(unique(significant_eCAVIAR_result_filter_snp_eQTL$gene_symble),
                                                        unique(significantly_fastENLOC_result_filter_snp_eQTL$gene_symble)))
    List_of_colocalizing_eGenes_intersection <- paste(intersect(unique(significant_eCAVIAR_result_filter_snp_eQTL$gene_symble),
                                                                unique(significantly_fastENLOC_result_filter_snp_eQTL$gene_symble)),collapse = ", ")

    significant_sGenes_intersection <- length(intersect(unique(significant_eCAVIAR_result_filter_snp_sQTL$gene_symble),
                                                        unique(significantly_fastENLOC_result_filter_snp_sQTL$gene_symble)))
    List_of_colocalizing_sGenes_intersection <- paste(intersect(unique(significant_eCAVIAR_result_filter_snp_sQTL$gene_symble),
                                                                unique(significantly_fastENLOC_result_filter_snp_sQTL$gene_symble)),collapse = ", ")

    significant_sGenes_by_intron_cluster_intersection <- length(intersect(unique(significant_eCAVIAR_result_filter_snp_sQTL$Gene_id_intron_cluster),
                                                                         unique(significantly_fastENLOC_result_filter_snp_sQTL$Gene_id_intron_cluster)))
    List_of_colocalizing_sGenes_by_intron_cluster_intersection <- paste(intersect(unique(significant_eCAVIAR_result_filter_snp_sQTL$Gene_id_intron_cluster),
                                                                               unique(significantly_fastENLOC_result_filter_snp_sQTL$Gene_id_intron_cluster)),
                                                                      collapse = ", ")

    significantly_colocalizing_eGenes_sGenes_intersection <- length(unique(c(intersect(unique(significant_eCAVIAR_result_filter_snp_eQTL$gene_symble),
                                                                                       unique(significantly_fastENLOC_result_filter_snp_eQTL$gene_symble)),
                                                                             intersect(unique(significant_eCAVIAR_result_filter_snp_sQTL$gene_symble),
                                                                                       unique(significantly_fastENLOC_result_filter_snp_sQTL$gene_symble)))))
    List_of_colocalizing_eGenes_sGenes_intersection <- paste(unique(c(intersect(unique(significant_eCAVIAR_result_filter_snp_eQTL$gene_symble),
                                                                                unique(significantly_fastENLOC_result_filter_snp_eQTL$gene_symble)),
                                                                      intersect(unique(significant_eCAVIAR_result_filter_snp_sQTL$gene_symble),
                                                                                unique(significantly_fastENLOC_result_filter_snp_sQTL$gene_symble)))),
                                                            collapse = ", ")
    Nearest_Gene <- dat_snps[rsID==snp,nearestGene]

    base_data <- eCAVIAR_result_filter_snp[1,c(1:11)]
    if (all(is.na(base_data))) {
      base_data <- subset(all_fastENLOC_result, rsID_GWAS == snp)[1,c(1:11)]
    }
    if (all(is.na(base_data))) {
      next
    }

    base_data$eQTL_tissue_geng_pairs_tested_eCAVIAR <- eQTL_tissue_geng_pairs_tested_eCAVIAR
    base_data$eGenes_tested_eCAVIAR <- eGenes_tested_eCAVIAR
    base_data$significant_eQTL_tissue_geng_pairs_tested_eCAVIAR <- significant_eQTL_tissue_geng_pairs_tested_eCAVIAR
    base_data$significant_eGenes_eCAVIAR <- significant_eGenes_eCAVIAR
    base_data$List_of_colocalizing_eGenes_eCAVIAR <- ifelse(is.null(List_of_colocalizing_eGenes_eCAVIAR),NA,List_of_colocalizing_eGenes_eCAVIAR)
    base_data$sQTL_tissue_geng_pairs_tested_eCAVIAR <- sQTL_tissue_geng_pairs_tested_eCAVIAR
    base_data$sGenes_tested_eCAVIAR <- sGenes_tested_eCAVIAR
    base_data$significant_sQTL_tissue_geng_pairs_tested_eCAVIAR <- significant_sQTL_tissue_geng_pairs_tested_eCAVIAR
    base_data$significant_sGenes_eCAVIAR <- significant_sGenes_eCAVIAR
    base_data$List_of_colocalizing_sGenes_eCAVIAR <- ifelse(is.null(List_of_colocalizing_sGenes_eCAVIAR),NA,List_of_colocalizing_sGenes_eCAVIAR)
    base_data$significant_sGenes_by_intron_cluster_eCAVIAR <- significant_sGenes_by_intron_cluster_eCAVIAR
    base_data$List_of_colocalizing_sGenes_by_intron_cluster_eCAVIAR <- ifelse(is.null(List_of_colocalizing_sGenes_by_intron_cluster_eCAVIAR),NA,List_of_colocalizing_sGenes_by_intron_cluster_eCAVIAR)
    base_data$eGenes_sGenes_tested_eCAVIAR <- eGenes_sGenes_tested_eCAVIAR
    base_data$significantly_colocalizing_eGenes_sGenes_eCAVIAR <- significantly_colocalizing_eGenes_sGenes_eCAVIAR
    base_data$List_of_significantly_colocalizing_eGenes_sGenes_eCAVIAR <- ifelse(is.null(List_of_significantly_colocalizing_eGenes_sGenes_eCAVIAR),NA,List_of_significantly_colocalizing_eGenes_sGenes_eCAVIAR)
    base_data$List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_eCAVIAR <- ifelse(is.null(List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_eCAVIAR),NA,List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_eCAVIAR)
    base_data$significantly_eQTL_tissue_geng_pairs_tested_fastENLOC <- significantly_eQTL_tissue_geng_pairs_tested_fastENLOC
    base_data$significantly_eGenes_tested_fastENLOC <- significantly_eGenes_tested_fastENLOC
    base_data$List_of_significantly_colocalizing_eGenes_fastENLOC <- ifelse(is.null(List_of_significantly_colocalizing_eGenes_fastENLOC),NA,List_of_significantly_colocalizing_eGenes_fastENLOC)
    base_data$significantly_sQTL_tissue_geng_pairs_tested_fastENLOC <- significantly_sQTL_tissue_geng_pairs_tested_fastENLOC
    base_data$significantly_sGenes_tested_fastENLOC <- significantly_sGenes_tested_fastENLOC
    base_data$List_of_significantly_colocalizing_sGenes_fastENLOC <- ifelse(is.null(List_of_significantly_colocalizing_sGenes_fastENLOC),NA,List_of_significantly_colocalizing_sGenes_fastENLOC)
    base_data$significant_sGenes_by_intron_cluster_fastENLOC <- significant_sGenes_by_intron_cluster_fastENLOC
    base_data$List_of_significantly_colocalizing_sGenes_by_intron_cluster_fastENLOC <- ifelse(is.null(List_of_significantly_colocalizing_sGenes_by_intron_cluster_fastENLOC),NA,List_of_significantly_colocalizing_sGenes_by_intron_cluster_fastENLOC)
    base_data$significantly_colocalizing_eGenes_sGenes_fastENLOC <- significantly_colocalizing_eGenes_sGenes_fastENLOC
    base_data$List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_fastENLOC <- ifelse(is.null(List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_fastENLOC),NA,List_of_tissues_with_significantly_colocalizing_esGenes_sGenes_fastENLOC)
    base_data$significant_eGenes_tested_union <- significant_eGenes_tested_union
    base_data$List_of_significant_colocalizing_eGenes_union <- ifelse(is.null(List_of_significant_colocalizing_eGenes_union),NA,List_of_significant_colocalizing_eGenes_union)
    base_data$significant_sGenes_tested_union <- significant_sGenes_tested_union
    base_data$List_of_significant_colocalizing_sGenes_union <- ifelse(is.null(List_of_significant_colocalizing_sGenes_union),NA,List_of_significant_colocalizing_sGenes_union)
    base_data$List_of_significant_colocalizing_sGenes_by_intron_cluster_union <- ifelse(is.null(List_of_significant_colocalizing_sGenes_by_intron_cluster_union),NA,List_of_significant_colocalizing_sGenes_by_intron_cluster_union)
    base_data$significant_sGenes_by_intron_cluster_union <- significant_sGenes_by_intron_cluster_union
    base_data$significantly_colocalizing_eGenes_sGenes_union <- significantly_colocalizing_eGenes_sGenes_union
    base_data$List_of_significantly_colocalizing_eGenes_sGenes_union <- ifelse(is.null(List_of_significantly_colocalizing_eGenes_sGenes_union),NA,List_of_significantly_colocalizing_eGenes_sGenes_union)
    base_data$List_of_significantly_colocalizing_tissues_from_eGenes_sGenes_union <- ifelse(is.null(List_of_significantly_colocalizing_tissues_from_eGenes_sGenes_union),NA,List_of_significantly_colocalizing_tissues_from_eGenes_sGenes_union)
    base_data$significant_eGenes_intersection <- significant_eGenes_intersection
    base_data$List_of_colocalizing_eGenes_intersection <- ifelse(is.null(List_of_colocalizing_eGenes_intersection),NA,List_of_colocalizing_eGenes_intersection)
    base_data$significant_sGenes_intersection <- significant_sGenes_intersection
    base_data$List_of_colocalizing_sGenes_intersection <- ifelse(is.null(List_of_colocalizing_sGenes_intersection),NA,List_of_colocalizing_sGenes_intersection)
    base_data$significant_sGenes_by_intron_cluster_intersection <- significant_sGenes_by_intron_cluster_intersection
    base_data$List_of_colocalizing_sGenes_by_intron_cluster_intersection <- ifelse(is.null(List_of_colocalizing_sGenes_by_intron_cluster_intersection),NA,List_of_colocalizing_sGenes_by_intron_cluster_intersection)
    base_data$significantly_colocalizing_eGenes_sGenes_intersection <- significantly_colocalizing_eGenes_sGenes_intersection
    base_data$List_of_colocalizing_eGenes_sGenes_intersection <- ifelse(is.null(List_of_colocalizing_eGenes_sGenes_intersection),NA,List_of_colocalizing_eGenes_sGenes_intersection)
    base_data$Nearest_Gene <- Nearest_Gene
    all_base_data <- rbind(base_data,all_base_data)

  }
  all_base_data <- all_base_data %>%  dplyr::mutate(across(where(is.character), ~ dplyr::na_if(., "")))

  fwrite(all_base_data,
         file=paste0("stable_summary_for_coloc_",trait,".txt"),
         na = "NA",
         quote = FALSE,
         sep = "\t")
}

sessionInfo()
