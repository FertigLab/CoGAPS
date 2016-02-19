#include "TransTest.h"

//[[Rcpp::export]]
Rcpp::List test_run(Rcpp::NumericVector y, Rcpp::IntegerVector treatStatus,
                    Rcpp::NumericVector timeRecorded, int iter) {
    TransTest obj(y, treatStatus, timeRecorded);

    obj.runMCMC(iter);

    return Rcpp::List::create(Rcpp::Named("beta0")=obj._beta0,
                              Rcpp::Named("beta1")=obj._beta1,
                              Rcpp::Named("tau")=obj._tau);
}
