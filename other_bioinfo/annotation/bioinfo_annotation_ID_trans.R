####################
# AIM:To tansfer your ID of bioinfo data
# AUTHOR:Cao Jun
# DATE:2023-03-24
###################

# load your packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(biomaRt)
library(AnnotationDbi)
library(AnnotationHub)
library(tidyverse)



# method 1 :AnnotationDbi

## choice1:online
ah <- AnnotationHub(localHub = FALSE)

## choice2:local
symbol <- mapIds(org.Mm.eg.db,
    keys = ensmel,
    column = "SYMBOL",
    keytype = "ENSEMBL", #' ENTREZID'
    multiVals = "first"
)


# method 2 :biomaRt
# ensembl_mart <- useEnsembl(biom)

## ensembl_transcript_id → external_gene_name
mart <- useMart("ensembl", "mmusculus_gene_ensembl")
## 人类选择hsapiens_gene_ensembl
gene_name <- getBM(
    attributes = c("ensembl_transcript_id", "external_gene_name", "ensembl_gene_id"),
    filters = "ensembl_transcript_id",
    values = gene,
    mart = mart
)



# method 3 :clusterProfiler::bitr
