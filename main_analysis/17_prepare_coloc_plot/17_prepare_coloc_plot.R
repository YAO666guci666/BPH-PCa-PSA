rm(list=ls())

gc()

getwd()

.libPaths()

library(data.table)

inputFile=c("ieu-b-4809_PC","finngen_BPH","GCST90476277_PSA")
width=10
height=8

all_Num_GWAS_loci_coloc <- data.frame()
for (trait in inputFile) {

  stable6 <- fread(paste0("../07_fuma_gwas/stable6.txt"))

  stable6_filter <- subset(stable6, Trait == trait)

  data <- fread(paste0("../16_summarize_coloc/stable_summary_for_coloc_",trait,".txt"))

  Num_GWAS_loci_coloc <- data.frame(Trait = trait)
  Num_GWAS_loci_coloc$Num_GWAS_loci <- length(unique(stable6_filter$rsID))

  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_eCAVIAR <- length(unique(subset(data,significant_eGenes_eCAVIAR>=1)$rsID_GWAS))
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_eCAVIAR_percentage <- Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_eCAVIAR/Num_GWAS_loci_coloc$Num_GWAS_loci
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_fastenloc <- length(unique(subset(data,significantly_eGenes_tested_fastENLOC>=1)$rsID_GWAS))
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_fastenloc_percentage <- Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_fastenloc/Num_GWAS_loci_coloc$Num_GWAS_loci

  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_eCAVIAR <- length(unique(subset(data,significant_sGenes_eCAVIAR>=1)$rsID_GWAS))
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_eCAVIAR_percentage <- Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_eCAVIAR/Num_GWAS_loci_coloc$Num_GWAS_loci
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_fastenloc <- length(unique(subset(data,significantly_sGenes_tested_fastENLOC>=1)$rsID_GWAS))
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_fastenloc_percentage <- Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_fastenloc/Num_GWAS_loci_coloc$Num_GWAS_loci

  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR <- length(unique(subset(data,significantly_colocalizing_eGenes_sGenes_eCAVIAR>=1)$rsID_GWAS))
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR_percentage <- Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR/Num_GWAS_loci_coloc$Num_GWAS_loci
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc <- length(unique(subset(data,significantly_colocalizing_eGenes_sGenes_fastENLOC>=1)$rsID_GWAS))
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc_percentage <- Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc/Num_GWAS_loci_coloc$Num_GWAS_loci

  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_union <- length(unique(subset(data,significant_eGenes_tested_union>=1)$rsID_GWAS))
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_union_percentage <- Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_union/Num_GWAS_loci_coloc$Num_GWAS_loci
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_union <- length(unique(subset(data,significant_sGenes_tested_union>=1)$rsID_GWAS))
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_union_percentage <- Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_union/Num_GWAS_loci_coloc$Num_GWAS_loci
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_union <- length(unique(subset(data,significantly_colocalizing_eGenes_sGenes_union>=1)$rsID_GWAS))
  Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_union_percentage <- Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_union/Num_GWAS_loci_coloc$Num_GWAS_loci

  all_Num_GWAS_loci_coloc <- rbind(all_Num_GWAS_loci_coloc, Num_GWAS_loci_coloc)

}

if (length(inputFile)>1) {

  Total_loci <- data.frame()

  Total_loci <- data.frame(Trait = "Total loci")
  Total_loci$Num_GWAS_loci <- sum(all_Num_GWAS_loci_coloc$Num_GWAS_loci)

  Total_loci$loci_greater_than_1_colocalizing_eGene_eCAVIAR <- sum(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_eCAVIAR)
  Total_loci$loci_greater_than_1_colocalizing_eGene_eCAVIAR_percentage <- Total_loci$loci_greater_than_1_colocalizing_eGene_eCAVIAR/Total_loci$Num_GWAS_loci
  Total_loci$loci_greater_than_1_colocalizing_eGene_fastenloc <- sum(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_fastenloc)
  Total_loci$loci_greater_than_1_colocalizing_eGene_fastenloc_percentage <- Total_loci$loci_greater_than_1_colocalizing_eGene_fastenloc/Total_loci$Num_GWAS_loci

  Total_loci$loci_greater_than_1_colocalizing_sGene_eCAVIAR <- sum(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_eCAVIAR)
  Total_loci$loci_greater_than_1_colocalizing_sGene_eCAVIAR_percentage <- Total_loci$loci_greater_than_1_colocalizing_sGene_eCAVIAR/Total_loci$Num_GWAS_loci
  Total_loci$loci_greater_than_1_colocalizing_sGene_fastenloc <- sum(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_fastenloc)
  Total_loci$loci_greater_than_1_colocalizing_sGene_fastenloc_percentage <- Total_loci$loci_greater_than_1_colocalizing_sGene_fastenloc/Total_loci$Num_GWAS_loci

  Total_loci$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR <- sum(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR)
  Total_loci$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR_percentage <- Total_loci$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR/Total_loci$Num_GWAS_loci
  Total_loci$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc <- sum(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc)
  Total_loci$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc_percentage <- Total_loci$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc/Total_loci$Num_GWAS_loci

  Total_loci$loci_greater_than_1_colocalizing_eGene_union <- sum(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_union)
  Total_loci$loci_greater_than_1_colocalizing_eGene_union_percentage <- Total_loci$loci_greater_than_1_colocalizing_eGene_union/Total_loci$Num_GWAS_loci
  Total_loci$loci_greater_than_1_colocalizing_sGene_union <- sum(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_union)
  Total_loci$loci_greater_than_1_colocalizing_sGene_union_percentage <- Total_loci$loci_greater_than_1_colocalizing_sGene_union/Total_loci$Num_GWAS_loci
  Total_loci$loci_greater_than_1_colocalizing_eGene_sGene_union <- sum(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_union)
  Total_loci$loci_greater_than_1_colocalizing_eGene_sGene_union_percentage <- Total_loci$loci_greater_than_1_colocalizing_eGene_sGene_union/Total_loci$Num_GWAS_loci

  Average <- data.frame()
  Average <- data.frame(Trait = "Average")
  Average$Num_GWAS_loci <- mean(all_Num_GWAS_loci_coloc$Num_GWAS_loci)

  Average$loci_greater_than_1_colocalizing_eGene_eCAVIAR <- mean(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_eCAVIAR)
  Average$loci_greater_than_1_colocalizing_eGene_eCAVIAR_percentage <- Average$loci_greater_than_1_colocalizing_eGene_eCAVIAR/Average$Num_GWAS_loci
  Average$loci_greater_than_1_colocalizing_eGene_fastenloc <- mean(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_fastenloc)
  Average$loci_greater_than_1_colocalizing_eGene_fastenloc_percentage <- Average$loci_greater_than_1_colocalizing_eGene_fastenloc/Average$Num_GWAS_loci

  Average$loci_greater_than_1_colocalizing_sGene_eCAVIAR <- mean(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_eCAVIAR)
  Average$loci_greater_than_1_colocalizing_sGene_eCAVIAR_percentage <- Average$loci_greater_than_1_colocalizing_sGene_eCAVIAR/Average$Num_GWAS_loci
  Average$loci_greater_than_1_colocalizing_sGene_fastenloc <- mean(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_fastenloc)
  Average$loci_greater_than_1_colocalizing_sGene_fastenloc_percentage <- Average$loci_greater_than_1_colocalizing_sGene_fastenloc/Average$Num_GWAS_loci

  Average$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR <- mean(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR)
  Average$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR_percentage <- Average$loci_greater_than_1_colocalizing_eGene_sGene_eCAVIAR/Average$Num_GWAS_loci
  Average$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc <- mean(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc)
  Average$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc_percentage <- Average$loci_greater_than_1_colocalizing_eGene_sGene_fastenloc/Average$Num_GWAS_loci

  Average$loci_greater_than_1_colocalizing_eGene_union <- mean(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_union)
  Average$loci_greater_than_1_colocalizing_eGene_union_percentage <- Average$loci_greater_than_1_colocalizing_eGene_union/Average$Num_GWAS_loci
  Average$loci_greater_than_1_colocalizing_sGene_union <- mean(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_sGene_union)
  Average$loci_greater_than_1_colocalizing_sGene_union_percentage <- Average$loci_greater_than_1_colocalizing_sGene_union/Average$Num_GWAS_loci
  Average$loci_greater_than_1_colocalizing_eGene_sGene_union <- mean(all_Num_GWAS_loci_coloc$loci_greater_than_1_colocalizing_eGene_sGene_union)
  Average$loci_greater_than_1_colocalizing_eGene_sGene_union_percentage <- Average$loci_greater_than_1_colocalizing_eGene_sGene_union/Average$Num_GWAS_loci

  all_Num_GWAS_loci_coloc <- rbind(all_Num_GWAS_loci_coloc, Average)
  all_Num_GWAS_loci_coloc <- rbind(all_Num_GWAS_loci_coloc, Total_loci)
}

fwrite(all_Num_GWAS_loci_coloc,
       file="all_Num_GWAS_loci_coloc.txt",
       sep="\t")

all_NumQTL_coloc_locus <- data.frame()
for (trait in inputFile) {

  stable6 <- fread(paste0("../07_fuma_gwas/stable6.txt"))

  stable6_filter <- subset(stable6, Trait == trait)

  data <- fread(paste0("../16_summarize_coloc/stable_summary_for_coloc_",trait,".txt"))

  trait_NumQTL_coloc_locus <- data.frame()

  for (snp in unique(data$rsID_GWAS)) {

  NumQTL_coloc_locus<- data.frame(leadSNPs = snp)

  NumQTL_coloc_locus$eGenes <- subset(data,rsID_GWAS==snp)$significant_eGenes_tested_union
  NumQTL_coloc_locus$sGenes <- subset(data,rsID_GWAS==snp)$significant_sGenes_tested_union
  NumQTL_coloc_locus$eGene_sGene_union <- subset(data,rsID_GWAS==snp)$significantly_colocalizing_eGenes_sGenes_union

  trait_NumQTL_coloc_locus <- rbind(trait_NumQTL_coloc_locus, NumQTL_coloc_locus)
  }

  trait_NumQTL_coloc_locus$Trait <- trait

  all_NumQTL_coloc_locus <- rbind(all_NumQTL_coloc_locus, trait_NumQTL_coloc_locus)

}

fwrite(all_NumQTL_coloc_locus,
       file="all_NumQTL_coloc_locus.txt",
       sep="\t")

sessionInfo()
