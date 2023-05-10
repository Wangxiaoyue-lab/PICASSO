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
