#ifndef _TRANS_H
#define _TRANS_H

#include <Rcpp.h>
#include <cmath>

// logistic transformation (log(p / 1 - p))
Rcpp::NumericVector logit(Rcpp::NumericVector data) {
    Rcpp::NumericVector transformation(data.size());
    double adj = Rcpp::max(data) + 1e-12;
    Rcpp::NumericVector y;
    y = transformation(y / (Rcpp::max(y) + 1e-12));

    for (int i = 0; i < data.size(); ++i) {
        transformation[i] = log(data[i] / (1 - data[i]));
    }

    return transformation;
}

// identity transformation (pattern is already linear)
Rcpp::NumericVector identity(Rcpp::NumericVector data) {
    return data;
}

#endif
