rm(list=ls())

gc()

getwd()

.libPaths()

library(Fast2TWAS)

significant_traits_file <- "../23_prepare_eclipser/demo_traits.txt"

ECLIPSER_py_install(is_help = T,
                    force = FALSE)

ECLIPSER_clumping_run(proxy_table = "../reference/ECLIPSER/GWASvar2gene_open_targets_proxy_table.tsv",
                      traits_table = "../reference/ECLIPSER/GWASvar2gene_open_targets_trait_table.tsv",
                      eQTL_table = "../reference/ECLIPSER/GWASvar2gene_open_targets_eQTL_table.tsv",
                      sQTL_table = "../reference/ECLIPSER/GWASvar2gene_open_targets_sQTL_table.tsv",
                       opentarget = "../reference/ECLIPSER/OpenTargetsGenetics_variant_and_gene_table_071620.csv",
                      significant_traits = significant_traits_file,

                      treat_sig_traits_as_one = NULL
                      )

sessionInfo()
