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
    double _mu0;    // beta prior mean
    double _tau0;   // beta prior precision
    double _a;      // variance prior shape
    double _b;      // variance prior rate

    // number of treatments
    int _nFactor;

    TransTest(arma::vec y, arma::ivec treatStatus, 
              arma::vec timeRecorded);

    ~TransTest() {};

    void update_pattern(std::vector<double>(*transformation)(std::vector<double>));

    void runMCMC(int iter);
};

#endif
