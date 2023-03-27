.libPaths("/public/home/luoliheng/miniconda3/envs/r1/lib/R/library")

# 使用monocle2进行分析
require(Seurat)
require(monocle)
require(tidyverse)
work_dir <- "/public/home/luoliheng/embryonic_cells/new_analysis/"
# Load the data
alldata <- readRDS(paste0(work_dir, "basic/qc.rds"))

# /public/home/luoliheng/embryonic_cells/new_analysis/trajectory/monocle2
cds <- as.CellDataSet(alldata) %>%
  detectGenes(min_expr = 0.1) %>%
  estimateSizeFactors() %>%
  estimateDispersions()

pdf(paste0(work_dir, "trajectory/monocle2/monocle2.pdf"))


cds <- setOrderingFilter(cds, VariableFeatures(alldata))
cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree", residualModelFormulaStr = "~type")


source("/public/home/luoliheng/embryonic_cells/new_analysis/trajectory/monocle2/order_cells.R")
require(igraph)
cds <- orderCells(cds)

saveRDS(cds, paste0(work_dir, "trajectory/monocle2/monocle2.rds"))

plot_ordering_genes(cds)
plot_cell_trajectory(cds, color_by = "type")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()


pseudotime_de <- differentialGeneTest(test[VariableFeatures(alldata), ],
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
write.csv(pseudotime_de, file = paste0(work_dir, "pseudotime_de.rds"), quote = FALSE, row.names = FALSE, col.names = TRUE)


states_de <- differentialGeneTest(test[expressed_genes, ],
  fullModelFormulaStr = "~State"
)
states_de <- states_de[order(states_de$qval), ]
write.csv(states_de, file = paste0(work_dir, "states_de.rds"), quote = FALSE, row.names = FALSE, col.names = TRUE)
