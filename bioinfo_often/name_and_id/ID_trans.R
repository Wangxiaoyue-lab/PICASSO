####################
#AIM:To tansfer your ID of bioinfo data
#AUTHOR:Cao Jun
#DATE:2023-03-24
###################

# load your packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(biomaRt)
library(AnnotationDbi)
library(AnnotationHub)
library(tidyverse)
library(conflicted)
conflict_scout()
conflict_prefer('filter','dplyr')

# method 1 :AnnotationDbi

## choice1:online 
ah <- AnnotationHub(localHub=FALSE)

## choice2:local
symbol <- mapIds(org.Mm.eg.db,
                keys=ensmel,
                column='SYMBOL',
                keytype='ENSEMBL',#'ENTREZID'
                multiVals='first')


# method 2 :biomaRt
ensembl_mart <- useEnsembl(biom)


# method 3 :clusterProfiler::bitr



