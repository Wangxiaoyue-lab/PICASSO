project_name <- "Embryo"
pal <- viridis::viridis(n = 10)
source("/public/home/luoliheng/SINGLE/R/utils.R", chdir = TRUE)
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(velocyto.R))
suppressPackageStartupMessages(require(SeuratWrappers))

work_dir <- "/public/home/luoliheng/embryonic_cells/new_analysis/rna_velo/"
project_name

velo_data <- ReadVelocity(file = "/public/home/luoliheng/embryonic_cells/rawdata/velocyto/combined.loom") %>%
  as.Seurat() %>%
  SCTransform(assay = "spliced") %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:20) %>%
  RunVelocity(deltaT = 1, kCells = 25, fit.quantile = 0.02)


# Add sample information to metadata
velo_data$type <- colnames(velo_data) %>%
  stringr::str_extract("^[^:]+")


ident.colors <- (scales::hue_pal())(n = length(x = levels(x = velo_data)))
names(x = ident.colors) <- levels(x = velo_data)
cell.colors <- ident.colors[Idents(object = velo_data)]
names(x = cell.colors) <- colnames(x = velo_data)

save_file(fun = saveRDS, name_string = "velo")
save_file(fun = pdf, name_string = "velo")
show.velocity.on.embedding.cor(
  emb = Embeddings(object = velo_data, reduction = "umap"), vel = Tool(
    object = velo_data,
    slot = "RunVelocity"
  ), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
  cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
  do.par = FALSE, cell.border.alpha = 0.1
)
dev.off()
