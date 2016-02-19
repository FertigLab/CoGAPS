#ifndef _TRANS_H
#define _TRANS_H

#include <Rcpp.h>
#include <cmath>

Rcpp::NumericVector logit(Rcpp::NumericVector data) {
    Rcpp::NumericVector transformation(data.size());

    for (int i = 0; i < data.size(); ++i) {
        transformation[i] = log(data[i] / (1 - data[i]));
    }

    return transformation;
}

#endif
