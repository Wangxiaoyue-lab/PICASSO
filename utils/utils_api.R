utils_wilcox_cpp <- function(x, y) {
    library(Rcpp)
    sourceCpp("./utils/cpp/fast_wilcox.cpp")
    wilcoxon_test(x, y)
}


utils_pca_cpp <- function(data, dims) {
    library(Rcpp)
    library(RcppArmadillo)
    sourceCpp("./utils/cpp/fast_pca.cpp")
    pca_cpp(data, dims)
}
