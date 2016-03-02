#ifndef _TRANSTEST_H
#define _TRANSTEST_H

#include <vector>
#include <cmath>
#include <Rcpp.h>

class TransTest {
  public:
    // stored data
    Rcpp::NumericVector y_all;
    Rcpp::IntegerVector _treatStatus;
    Rcpp::NumericVector _timeRecorded;

    // current estimates
    Rcpp::NumericVector _beta0;
    Rcpp::NumericVector _beta1;
    Rcpp::NumericVector _tau;

    // chains
    Rcpp::NumericMatrix _beta1chain;

    // priors for Gibbs Sampling of regression coefficients
    double _mu0;    // beta prior mean
    double _tau0;   // beta prior precision
    double _a;      // variance prior shape
    double _b;      // variance prior rate

    // number of treatments
    int _nFactor;

    TransTest(Rcpp::NumericVector y, Rcpp::IntegerVector treatStatus, 
              Rcpp::NumericVector timeRecorded);

    ~TransTest() {};

    void update_pattern(Rcpp::NumericVector(*transformation)(Rcpp::NumericVector));

    void runMCMC(int iter);
};

#endif
