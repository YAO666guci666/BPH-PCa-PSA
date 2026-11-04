library(data.table)
library(snpStats)
library(LAVA)
library(parallel)

loci = read.loci("./support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")

pairs_to_run = list(c('PSA','BPH'),c('PSA','PC'))
for (traits in pairs_to_run) {
input = process.input(
  input.info.file = "input.info.txt",
  sample.overlap.file = NULL,
  ref.prefix = "EUR",
  phenos = traits
)
trait_pair = paste(traits, collapse = "_")

main_output_dir = paste0("./", trait_pair)
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
  message("Creating the main output directory:", main_output_dir)
}

for (trait in traits) {
  trait_dir = file.path(main_output_dir, trait)
  if (!dir.exists(trait_dir)) {
    dir.create(trait_dir)
    message("Creating directory:", trait_dir)
  }
}

pair_dir = file.path(main_output_dir, trait_pair)
if (!dir.exists(pair_dir)) {
  dir.create(pair_dir)
  message("Creating directory:", pair_dir)
}

for (i in 1:nrow(loci)) {
  locus = process.locus(loci[i,], input)
  if (is.null(locus)) {
    next
  } else {
    tmp <- run.univ.bivar(locus, univ.thresh = 0.05/nrow(loci))

    result <- as.data.frame(tmp$univ)
    result <- cbind(result, data.frame(LOC = i))

    for (j in 1:length(traits)) {
      trait_data <- result[j, ]
      trait <- traits[j]
      trait_dir = file.path(main_output_dir, trait)

      save(trait_data, file = file.path(trait_dir, paste0(trait, "_loc", i, ".RData")))

    }

    bivar_data <- tmp$bivar
    if (is.null(bivar_data)) {

      bivar_data <- data.frame()
    } else {
      bivar_data <- cbind(as.data.frame(bivar_data), data.frame(LOC = i))
    }

    save(bivar_data, file = file.path(pair_dir, paste0(trait_pair, "_loc", i, ".RData")))

    if (i %% 100 == 0) {
      message("Processed locus", i, "/", nrow(loci), " (", round(i/nrow(loci)*100, 2), "%)")
    }
  }
}

merge_results <- function(trait, main_dir) {
  trait_dir = file.path(main_dir, trait)
  files = list.files(trait_dir, full.names = TRUE, pattern = '*.RData')

  all_data = NULL
  for (file in files) {
    load(file)

    if (exists("trait_data")) {
      data = trait_data

      if (is.null(all_data)) {
        all_data = data
      } else {
        all_data = rbind(all_data, data)
      }
    }
  }

  if (!is.null(all_data) && nrow(all_data) > 0) {
    all_data = merge(all_data, loci, by = 'LOC')
  }

  return(all_data)
}

for (trait in traits) {
  cat("Merging", trait, "results...\n")
  trait_results = merge_results(trait, main_output_dir)

  if (!is.null(trait_results) && nrow(trait_results) > 0) {
    output_file = file.path(main_output_dir, trait, paste0(trait, "_all_result.txt"))
    write.table(trait_results, file = output_file, sep = '\t', row.names = FALSE)
    cat("Saved", trait, "results Info:", output_file, "\n")
  } else {
    cat("No", trait, "valid results \n")
  }
}

cat("Univariate results merged!\n")

merge_bivar_results <- function(trait_pair, main_dir) {
  pair_dir = file.path(main_dir, trait_pair)
  files = list.files(pair_dir, full.names = TRUE, pattern = '*.RData')

  all_data = NULL
  for (file in files) {
    load(file)

    if (exists("bivar_data") && !is.null(bivar_data) && nrow(bivar_data) > 0) {
      if (is.null(all_data)) {
        all_data = bivar_data
      } else {
        all_data = rbind(all_data, bivar_data)
      }
    }
  }

  if (!is.null(all_data) && nrow(all_data) > 0) {
    all_data = merge(all_data, loci, by = 'LOC')
  }

  return(all_data)
}

cat("Merging bivariate results...\n")
bivar_results = merge_bivar_results(trait_pair, main_output_dir)

if (!is.null(bivar_results) && nrow(bivar_results) > 0) {
  output_file = file.path(main_output_dir, paste0(trait_pair, "_all_result.txt"))
  write.table(bivar_results, file = output_file, sep = '\t', row.names = FALSE)
  cat("Saved Info results Info:", output_file, "\n")

  if ("p" %in% colnames(bivar_results)) {
    cat("Computing FDR Info Filtering significant results...\n")
    fdr = p.adjust(bivar_results$p, method = "BH")
    fdr_results = data.frame(bivar_results, FDR = fdr)
    significant_results = subset(fdr_results, FDR < 0.05)

    fdr_file = file.path(main_output_dir, trait_pair, paste0(trait_pair, "_all_result_fdr.txt"))
    write.table(significant_results, file = fdr_file, sep = '\t', row.names = FALSE)
    cat("Saved FDR Info results (", nrow(significant_results), "significant loci) Info:", fdr_file, "\n")
  } else {
    warning("No 'p' columns, Cannot calculate FDR")
  }
} else {
  cat("No Info valid results \n")
}

cat("Bivariate results merged!\n")
}
sessionInfo()
