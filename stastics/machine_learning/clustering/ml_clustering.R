# 默认方向为--列为样本，行为特征
# 层次聚类及剪枝
ml_cluster_hclust <- function(object, ncores = NULL, cut = F, k = NULL, ...) {
    ncores <- ncores %||% 20
    hcc <- hclust(
        parallelDist::parDist(t(object),
            threads = ncores,
            method = "euclidean"
        ),
        method = "ward.D2"
    )
    if (cut == F) {
        return(hcc)
    } else {
        hc.cluster <- cutree(hcc, k...)
    }
}



ml_cluster <- function(object, methods, ...) {
    ml_cluster_kmeans <- function(...) {
        next
    }
    ml_cluster_kmedoids <- function(...) {
        next
    }
    ml_cluster_kmods <- function(...) {
        next
    }
    ml_cluster_krototypes <- function(...) {
        next
    }
    ml_cluster_fcm <- function(...) {
        next
    }
    ml_cluster_pam <- function(...) {
        next
    }
    ml_cluster_dbscan <- function(...) {
        next
    }
    ml_cluster_hdbscan <- function(...) {
        next
    }
    ml_cluster_gmm <- function(...) {
        next
    }
    ml_cluster_spectral <- function(...) {
        next
    }
    ml_cluster_ap <- function(...) {
        next
    }
}

# 双聚类
ml_bicluster <- function(...) {
    next
}

# 一致性聚类
ml_ccp <- function(..) {
    next
}
