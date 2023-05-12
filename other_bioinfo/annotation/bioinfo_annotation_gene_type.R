library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v86)


# EnsDb.Hsapiens.v79
# EnsDb.Hsapiens.v75

# EnsDb.Hsapiens.v79
# EnsDb.Hsapiens.v75

edb <- EnsDb.Hsapiens.v86
keys <- keys(edb, keytype = "GENEID")
gene2sym <- select(edb,
    keys = keys,
    columns = c("SYMBOL", "ENTREZID", "GENEBIOTYPE", "GENENAME"),
    keytype = "GENEID"
)
table(gene2sym$GENEBIOTYPE)

### 编码RNA
mRNA <- "protein_coding"

### 非编码RNA[最可能是lncRNA的，不是以基因定位来选择]
ncRNA <- c("sense_overlapping", "lincRNA", "3prime_overlapping_ncRNA", "processed_transcript", "sense_intronic", "bidirectional_promoter_lncRNA", "non_coding", "antisense_RNA")
