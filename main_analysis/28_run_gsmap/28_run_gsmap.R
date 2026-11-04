rm(list=ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)

inputFile=c("GCST90476277_PSA","finngen_BPH","ieu-b-4809_PC")

hdf5_path="../27_prepare_gsmap/ST/E16.5_E1S1.MOSTA.h5ad"
annotation='annotation'
data_layer='count'

for (trait in inputFile) {
  gsMap_format_sumstats_run(sumstats=paste0("../02_harmonize_gwas/",trait,"_hg38_qc.txt"),
                            snp="SNP",
                            a1="effect_allele",
                            a2="other_allele",
                            info=NULL,
                            beta="beta",
                            se="se",
                            p="pval",
                            frq="eaf",
                            n="N",
                            chr="CHR",
                            pos="BP",
                            dbsnp=NULL,
                            chunksize="500000",
                            format="gsMap",
                            info_min=NULL,
                            maf_min=NULL,
                            outname=trait)

}

for (trait in inputFile) {
  gsMap_quick_mode_run(sample_name=gsub(".h5ad","",basename(hdf5_path)),
                       gsMap_resource_dir='gsMap_resource',
                       hdf5_path=hdf5_path,
                       annotation=annotation,
                       data_layer=data_layer,
                       trait_name=trait,
                       sumstats_file=paste0('./result/',trait,".sumstats.gz"),
                       sumstats_config_file=NULL,
                       homolog_file='./gsMap_resource/homologs/mouse_human_homologs.txt',
                       max_processes=parallel::detectCores()/2,
                       latent_representation=NULL,
                       num_neighbour=NULL,
                       num_neighbour_spatial=NULL,
                       gM_slices=NULL,
                       is_pearson_residuals=FALSE,
                       argsfile=NULL)
}

sessionInfo()
