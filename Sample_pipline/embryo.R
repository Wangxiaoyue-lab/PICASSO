suppressPackageStartupMessages(require(velocyto.R))


work_dir <- "/public/home/luoliheng/embryonic_cells/new_analysis/trajectory/rna_velo/"
project_name <- "Embryo"

source("/public/home/luoliheng/SINGLE/scRNA_analysis/utils.R", chdir = TRUE)

velo_data <- readRDS(save_file(fun = saveRDS, name_string = "velo_noharmony"))

save_file(fun = pdf, name_string = "velo_noharmony")

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = velo_data)))
names(x = ident.colors) <- levels(x = velo_data)
cell.colors <- ident.colors[Idents(object = velo_data)]
names(x = cell.colors) <- colnames(x = velo_data)


print(show.velocity.on.embedding.cor(
  emb = Embeddings(object = velo_data, reduction = "umap"), vel = Tool(
    object = velo_data,
    slot = "RunVelocity"
  ), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
  cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
  do.par = FALSE, cell.border.alpha = 0.1
))
dev.off()
