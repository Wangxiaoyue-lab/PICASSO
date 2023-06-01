sc_batch_kbet <- function(object, batch_label, data_select) {
    # https://github.com/theislab/kBET
    library(kBET)
    library(FNN)
    assertthat::assert_that(class(object) == "Seurat")
    assertthat::assert_that(batch_label %in% colnames(object@meta.data))
    if (data_select == "pca") {
        # data <- as.data.frame(object@reductions$pca@cell.embeddings)
    } else if (data_select == "harmony") {
        # data <- as.data.frame(object@dr$harmony@cell.embeddings)
    } else if (data_select == "cca") {

    }
    # data: a matrix (rows: samples, columns: features (genes))
    batch <- object@meta.data %>%
        pull(!!sym(batch_label))
    k0 <- floor(mean(table(batch))) # neighbourhood size: mean batch size
    knn <- get.knn(data, k = k0, algorithm = "cover_tree")
    batch.estimate <- kBET(data, batch, k = k0, knn = knn, plot = TRUE)
}

sc_batch_silhouette <- function() {
    # https://github.com/theislab/kBET
    batch.silhouette <- batch_sil(pca.data, batch)
}


sc_batch_pcRegression <- function() {
    # https://github.com/theislab/kBET
    batch.pca <- pcRegression(pca.data, batch)
}


sc_batch_lisi <- function() {
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6964114/
}
