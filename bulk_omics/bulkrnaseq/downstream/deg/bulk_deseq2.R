library(DESeq2)
library(tximport)

# 1 Accepted data
## The count data is rocommanded


## Other data
# Salmon
# Sailfish
# kallisto
# RSEM
# eXpress
# cufflinks/cuffdiff

txim <- tximport(files, type = "rsem", txOut = F)
# c("none", "salmon", "sailfish", "alevin", "kallisto", "rsem", "stringtie"




# 2 Create DDS object
dds <- DESeqDataSetFromTximport(txim, colData = samples, design = ~condition)


# 3 EDA



# 4 DEG
dds <- DESeq(dds)
res <- results(dds)
