rm(list=ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)

inputFile=c("GCST90476277_PSA","finngen_BPH")

hdf5_path="../27_prepare_gsmap/ST/E16.5_E1S1.MOSTA.h5ad"

gsMap_py_install(is_help = FALSE,
                 subcommand=NULL,
                 force = FALSE)

sessionInfo()
