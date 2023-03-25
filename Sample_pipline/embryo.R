require(Seurat)
require(monocle)
require(tidyverse)
work_dir <- "/public/home/luoliheng/embryonic_cells/new_analysis/"
# Load the data
alldata <- readRDS(paste0(work_dir, "basic/qc.rds"))

cds <- as.CellDataSet(alldata) %>%
  detectGenes(min_expr = 0.1) %>%
  estimateSizeFactors() %>%
  estimateDispersions()


cds <- setOrderingFilter(VariableFeatures(alldata)) %>%
  plot_ordering_genes() %>%
  reduceDimension(max_components = 2, method = "DDRTree") %>%
  orderCells()

pdf(paste0(work_dir, "monocle2/monocle2.pdf"))
plot_cell_trajectory(cds, color_by = "type")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()

saveRDS(cds, paste0(work_dir, "monocle2/monocle2.rds"))


pseudotime_de <- differentialGeneTest(test[VariableFeatures(alldata), ],
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]

states_de <- differentialGeneTest(test[expressed_genes, ],
  fullModelFormulaStr = "~State"
)
states_de <- states_de[order(states_de$qval), ]
