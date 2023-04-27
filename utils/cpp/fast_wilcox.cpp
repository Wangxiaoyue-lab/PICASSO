#include <Rcpp.h>
#include <vector>
#include <utility>
#include <algorithm>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double wilcoxon_test(NumericVector x, NumericVector y) {
    int n1 = x.size();
    int n2 = y.size();
    vector<pair<double, int>> all_values;
    for (int i = 0; i < n1; i++) {
        all_values.push_back(make_pair(x[i], 1));
    }
    for (int i = 0; i < n2; i++) {
        all_values.push_back(make_pair(y[i], 2));
    }
    sort(all_values.begin(), all_values.end());
    double rank_sum = 0;
    for (int i = 0; i < n1 + n2; i++) {
        if (all_values[i].second == 1) {
            rank_sum += i + 1;
        }
    }
    double U = rank_sum - n1 * (n1 + 1) / 2;
    double mean = n1 * n2 / 2;
    double sd = sqrt(n1 * n2 * (n1 + n2 + 1) / 12);
    return R::pnorm((U - mean) / sd,0.0,1.0,0,0);
}