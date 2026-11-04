library(data.table)

pairs <- c("BPH_PSA", "BPH_PCA", "PCA_PSA")

for(p in pairs){

  f_in <- file.path("../31_placo/PLACO_pairwise", p, "placo_dat.txt")
  if(!file.exists(f_in)) stop("Missing PLACO input:", f_in)

  d <- fread(f_in)
  setnames(d, tolower(names(d)))

  if("chr" %in% names(d)) setnames(d, "chr", "CHR")
  if("bp" %in% names(d)) setnames(d, "bp", "BP")
  if("pos" %in% names(d)) setnames(d, "pos", "BP")
  if("p.placo" %in% names(d)) setnames(d, "p.placo", "P")
  if("p" %in% names(d) && !"P" %in% names(d)) setnames(d, "p", "P")

  need <- c("CHR", "BP", "P")
  if(!all(need %in% names(d)))
    stop("placo_dat Missing required columns (", paste(setdiff(need, names(d)), collapse = ", "), "): ", f_in,
         "\n Actual columns:", paste(names(d), collapse = ", "))

  d <- d[, .(chromosome = CHR, base_pair_location = BP, p_value = P)]

  if(any(d$base_pair_location %% 1 != 0)) stop(p, "BP Info, check Info columns")

  n0 <- nrow(d)
  d <- d[chromosome %in% 1:22]
  if(nrow(d) < n0) cat(p, "Removed non- 1-22 chromosomes:", n0 - nrow(d), "rows \n")

  d <- d[!is.na(p_value) & p_value >= 0 & p_value <= 1]

  n_zero <- sum(d$p_value == 0)
  if(n_zero > 0){
    d$p_value[d$p_value == 0] <- min(d$p_value[d$p_value > 0], na.rm = TRUE)
    cat(p, "P=0 underflow handling:", n_zero, "rows \n")
  }

  setorder(d, p_value)
  n_dup <- sum(duplicated(d, by = c("chromosome", "base_pair_location")))
  d <- unique(d, by = c("chromosome", "base_pair_location"))
  if(n_dup > 0) cat(p, "Duplicated coordinates removed:", n_dup, "rows \n")

  f_out <- paste0("./FUMA_input_", p, ".txt")
  fwrite(d, f_out, sep = "\t")
  cat(p, "->", f_out, "| rows:", nrow(d),
      "| header: chromosome, base_pair_location, p_value\n")
}

cat("\n========================================\n")
cat("Submit FUMA Info Note:\n")
cat("1) Select 'Input file is in GRCh38'\n")
cat("2) Sample size (N) parameter: BPH_PSA=596163, BPH_PCA=377335, PCA_PSA=584078\n")
cat("3) Use the same reference panel 1000G EUR\n")
cat("========================================\n")

sessionInfo()
