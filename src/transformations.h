#ifndef _TRANS_H
#define _TRANS_H

#include <Rcpp.h>
#include <cmath>

// logistic transformation (log(p / 1 - p))
Rcpp::NumericVector logit(Rcpp::NumericVector data) {
    Rcpp::NumericVector transformation(Rcpp::clone(data));

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
