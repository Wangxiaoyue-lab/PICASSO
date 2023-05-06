library(immunedeconv)
# https://omnideconv.org/immunedeconv/

# quantiseq  10
# timer  6
# cibersort  22
# cibersort_abs
# mcp_counter    8
# xcell  64
# epic   6
# abis
# consensus_tme
# estimate

immunedeconv::deconvolute(gene_expression = gene_expression, methods = "")
