getwd()

rm(list=ls())

gc()

getwd()

.libPaths()

if (!dir.exists("result")) {
  dir.create("result")
}

if (!requireNamespace("Seurat", quietly = TRUE))install.packages("Seurat")
if (!requireNamespace("Signac", quietly = TRUE))install.packages("Signac")
if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE))BiocManager::install("EnsDb.Hsapiens.v75")
if (!requireNamespace("Open4Gene", quietly = TRUE))devtools::install_github("hbliu/Open4Gene")
if (!requireNamespace("bedtoolsr", quietly = TRUE))remotes::install_github("PhanstielLab/bedtoolsr")

library(Seurat)
library(Signac)
library(Fast2TWAS)
library(Open4Gene)
library(bedtoolsr)
library(EnsDb.Hsapiens.v75)

load("./result/Open4Gene_input.Rdata")

RNA.Counts=pbmc.rna@assays[["RNA"]]@layers[["counts"]]
ATAC.Counts=pbmc.atac@assays[["ATAC"]]@counts
rownames(RNA.Counts) <- rownames(pbmc.rna)
colnames(RNA.Counts) <- colnames(pbmc.rna)
Meta.Data <- pbmc.rna@meta.data

all_cols <- colnames(Meta.Data)
find_col <- function(candidates) {
  for (c in candidates) {
    if (c %in% all_cols) return(c)
  }
  return(NULL)
}
cov_rna <- find_col(c("log10_nCount_RNA", "nCount_RNA", "log10_nUMI", "nUMI"))
cov_mt  <- find_col(c("percent.mt", "percent.MT", "percent_mito", "mito_percent"))
Covariates <- c(cov_rna, cov_mt)
Covariates <- Covariates[!sapply(Covariates, is.null)]
if (length(Covariates) == 0) Covariates <- NULL
cat("Info Using Covariates:", paste(Covariates, collapse = ", "), "\n")

Celltype_col="seurat_annotations"
Peak2Gene.Dis=100000
Celltype="All"
Binary = FALSE
MinNum.Cells=1

Gene.Annotation <- genes(EnsDb.Hsapiens.v75)

std_chrs <- c(as.character(1:22), "X", "Y", "MT")
Gene.Annotation <- Gene.Annotation[
  as.character(seqnames(Gene.Annotation)) %in% std_chrs
]

ga_df <- as.data.frame(Gene.Annotation)
ga_df$width <- ga_df$end - ga_df$start + 1
ga_df <- ga_df[order(ga_df$gene_name, -ga_df$width), ]
ga_df <- ga_df[!duplicated(ga_df$gene_name), ]

Gene.Annotation <- GenomicRanges::makeGRangesFromDataFrame(
  ga_df, keep.extra.columns = TRUE
)

options(timeout = 10000)

bedtools.path <- Open4Gene_install(force=T)
options(bedtools.path =dirname(bedtools.path))

result <- Open4Gene_run(RNA = RNA.Counts,
                        ATAC = ATAC.Counts,
                        Meta.data = Meta.Data,
                        Covariates=Covariates,
                        Celltype_col=Celltype_col,

                        Peak2Gene.Dis=Peak2Gene.Dis,
                        Gene.Annotation = Gene.Annotation,
                        Res=NULL,
                        Celltype = Celltype,
                        Binary = Binary,
                        Method = "hurdle",
                        MinNum.Cells = MinNum.Cells)

write.table(result[["Open4Gene.result"]]@Res,
            file = "Open4Gene.obj.res.txt",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)

save(result, file = "Open4Gene.obj.Rdata")

sessionInfo()
