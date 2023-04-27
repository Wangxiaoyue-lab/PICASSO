#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List pca_cpp(const mat& data, int ncomp) {
    int n = data.n_rows, p = data.n_cols;
    mat X = data;
    rowvec colmeans = mean(X);
    X.each_row() -= colmeans;
    mat U, V;
    vec s;
    svd(U, s, V, X);
    mat Y = U.cols(0, ncomp - 1) * diagmat(s.subvec(0, ncomp - 1));
    return List::create(Named("scores") = Y,
                        Named("loadings") = V.cols(0, ncomp - 1),
                        Named("sdev") = s.subvec(0, ncomp - 1));
}