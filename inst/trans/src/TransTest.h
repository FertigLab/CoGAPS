#ifndef _TRANSTEST_H
#define _TRANSTEST_H

#include <vector>
#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

class TransTest {
  public:
    // stored data
    arma::vec y_all;
    arma::ivec _treatStatus;
    arma::vec _timeRecorded;

    // current estimates
    arma::vec _beta0;
    arma::vec _beta1;
    arma::vec _tau;

    // priors for Gibbs Sampling of regression coefficients
    double _mu0 = 0.0;    // beta prior mean
    double _tau0 = 0.001; // beta prior precision
    double _a = 1.0;      // variance prior shape
    double _b = 1.0;      // variance prior rate

    // number of treatments
    int _nFactor = 2;

    TransTest(arma::vec y, arma::ivec treatStatus, 
              arma::vec timeRecorded);

    ~TransTest() {};

    void update_pattern(std::vector<double>(*transformation)(std::vector<double>));

    void runMCMC(int iter);
};

#endif
