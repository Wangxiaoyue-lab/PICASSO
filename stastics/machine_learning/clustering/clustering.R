ml_cluster <- function(object, methods) {
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

# https://github.com/LuyiTian/sc_mixology/blob/master/script/clustering/cluster_eval.R 评估聚类结果

cal_entropy <- function(x) {
    x <- x / sum(x)
    -sum(x * log(x))
}

ARI_matric <- function(x, y) {
    x <- factor(x)
    y <- factor(y)
    return(adjustedRandIndex(x, y))
}

ml_cluster_metric <- function(..., methods) {
    ml_cluster_metric_ami <- function(...) {
        next
    }
    ml_cluster_metric_ari <- function(...) {
        next
    }
    ml_cluster_metric_psi <- function(...) {
        next
    }
    ml_cluster_metric_nmi <- function(...) {
        next
    }
}
