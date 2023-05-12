# 1 from transcript ID to gene ID
# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

library(tximport)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


# Salmon
# Sailfish
# kallisto
# RSEM
# eXpress
# cufflinks/cuffdiff


# 准备tx2gene文件
## 这是官方推荐方法，但是好像行不通了，提出来不对
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
write.csv(tx2gene, file = "/other_bioinfo/annotation/annotation_files/reference/tx2gene_hg38.csv")




# 2 from transcript matrix to gene matrix
# https://www.jianshu.com/p/e0acb957b351

files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
all(file.exists(files))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
