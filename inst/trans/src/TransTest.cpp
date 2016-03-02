#include "TransTest.h"
#include "transformation.h"

TransTest::TransTest(Rcpp::NumericVector y, 
                     Rcpp::IntegerVector treatStatus, 
                     Rcpp::NumericVector timeRecorded) {
    y_all = y;
    _treatStatus = treatStatus;
    _timeRecorded = timeRecorded;

    // priors for Gibbs Sampling of regression coefficients
    _mu0 = 0.0;    // beta prior mean
    _tau0 = 0.001; // beta prior precision
    _a = 1.0;      // variance prior shape
    _b = 1.0;      // variance prior rate


    // number of treatments
    _nFactor = 2;

    // initialize parameters from prior
    _beta0 = Rcpp::rnorm(_nFactor, _mu0, sqrt(1. / _tau0));
    _beta1 = Rcpp::rnorm(_nFactor, _mu0, sqrt(1. / _tau0));
    _tau = Rcpp::rgamma(_nFactor, _a, _b);
}

void TransTest::update_pattern(Rcpp::NumericVector(*transformation)(Rcpp::NumericVector)) {
    // http://www.cs.toronto.edu/~radford/csc2541.S11/week3.pdf

    for (int i = 0; i < _nFactor; ++i) {
        // initialize lists of y and x for each regression
        Rcpp::NumericVector y = y_all[_treatStatus == i];
        y = transformation(y / (Rcpp::max(y) + 1e-12));
        Rcpp::NumericVector x = _timeRecorded[_treatStatus == i];
        Rcpp::NumericVector x_sq = Rcpp::pow(x, 2.0);

        // initialize variables for full conditionals
        double post_mean, post_var;     // normal distributions
        double post_shape, post_rate;   // gamma distributions
        int n = y.size();

        // \beta_0 | x, y, \beta_1, \tau
        post_var = 1. / (_tau0 + n * _tau[i]);
        post_mean = (_tau0 * _mu0 + _tau[i] * Rcpp::sum(y - _beta1[i] * x)) * post_var;
        _beta0[i] = Rcpp::as<double>(Rcpp::rnorm(1, post_mean, 
                                                 sqrt(post_var)));
        // \beta_1 | x, y, \beta_0, \tau
        post_var = 1. / (_tau0 + _tau[i] * Rcpp::sum(x_sq));
        post_mean = (_tau0 * _mu0 + _tau[i] * Rcpp::sum(x * (y - _beta1[i]))) * post_var;
        _beta1[i] = Rcpp::as<double>(Rcpp::rnorm(1, post_mean, 
                                                 sqrt(post_var)));

        // \tau | x, y, _beta0, _beta1
        post_shape = _a + n / 2.;
        post_rate = _b + Rcpp::sum(Rcpp::pow(y - _beta0[i] - _beta1[i] * x, 2.0) / 2.0);
        _tau[i] = Rcpp::as<double>(Rcpp::rgamma(1, post_shape, 
                                                1. / post_rate));
    }
}

void TransTest::runMCMC(int iter) {
    // size chains
    Rcpp::NumericMatrix chain(iter, _nFactor);

    for (int s = 0; s < iter; ++s) {
        update_pattern(&logit);
        for (int j = 0; j < _nFactor; ++j) {
            chain(s, j) = _beta1[j];
        }
    }

    _beta1chain = chain;
}
