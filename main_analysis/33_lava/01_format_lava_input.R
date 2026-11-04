library(data.table)
library(tidyr)

input_dir <- "./input"
output_dir <- "./formatted"
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
files <- c(BPH="BPH37.txt",PC="PC37.txt",PSA="PSA37.txt")
for (pheno in names(files)) {
  file_path <- file.path(input_dir,files[pheno])
  if(!file.exists(file_path)) stop("Missing LAVA input:",file_path)
  data <- fread(file_path)
  print(names(data))
  bcol <- intersect(c("b","beta","BETA","Effect"),names(data))[1]
  secol <- intersect(c("se","SE","StdErr"),names(data))[1]
  ncol0 <- intersect(c("n","N"),names(data))[1]
  if(any(is.na(c(bcol,secol,ncol0)))) stop(pheno,"Missing beta/se/N columns")
  data$b <- as.numeric(data[[bcol]])
  data$se <- as.numeric(data[[secol]])
  data$Z <- data$b/data$se
  data$n_lava <- as.numeric(data[[ncol0]])
  f <- data[,c("SNP","A1","A2","Z","n_lava")]
  colnames(f) <- c("SNP","A1","A2","Z","N")
  fwrite(f,file.path(output_dir,paste0(pheno,".txt")),sep="	")
}
sessionInfo()
