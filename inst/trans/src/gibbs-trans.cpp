#include <vector>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

std::vector<double> logit(std::vector<double> data) {
    std::vector<double> transformation(data.size());

    for (int i = 0; i < data.size(); ++i) {
        transformation[i] = log(data[i] / (1 - data[i]));
    }

    return transformation;
}


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
              arma::vec timeRecorded) {
        y_all = y;
        _treatStatus = treatStatus;
        _timeRecorded = timeRecorded;

        _beta0.resize(2);
        _beta1.resize(2);
        _tau.resize(2);
    }

    void update_pattern(std::vector<double>(*transformation)(std::vector<double>)) {
        // http://www.cs.toronto.edu/~radford/csc2541.S11/week3.pdf

        for (int i = 0; i < _nFactor; ++i) {
            // initialize lists of y and x for each regression
            arma::vec y = y_all(arma::find(_treatStatus == i));
            arma::vec x = _timeRecorded(arma::find(_treatStatus == i));
            arma::vec x_sq(x.size());

            for (int j = 0; j < x.size(); ++j) {
                x_sq[j] = pow(x[j], 2.0);
            }

            // initialize variables for full conditionals
            double post_mean, post_var;     // normal distributions
            double post_shape, post_rate;   // gamma distributions
            int n = y.size();

            // \beta_0 | x, y, \beta_1, \tau
            post_var = 1. / (_tau0 + n * _tau[i]);
            post_mean = (_tau0 * _mu0 + _tau[i] * arma::sum(y - _beta1[i] * x)) * post_var;
            _beta0[i] = arma::conv_to<double>::from(arma::randn(1) * sqrt(post_var) + post_mean);

            // \beta_1 | x, y, \beta_0, \tau
            post_var = 1. / (_tau0 + _tau[i] * arma::sum(x_sq));
            post_mean = (_tau0 * _mu0 + _tau[i] * arma::sum(x * (y - _beta1[i]))) * post_var;
            _beta1[i] = arma::conv_to<double>::from(arma::randn(1) * sqrt(post_var) + post_mean);

            // \tau | x, y, _beta0, _beta1
            post_shape = _a + n / 2.;
            post_rate = _b + arma::sum(arma::pow(y - _beta0[i] - _beta1[i] * x, 2.0) / 2.0);
            _tau[i] = arma::conv_to<double>::from(arma::randg(1, arma::distr_param(post_shape, post_rate)));
        }
    }

    void runMCMC(int iter) {
        for (int s = 0; s < iter; ++s) {
            update_pattern(&logit);
        }
    }
};

//[Rcpp::export]]
Rcpp::List test_run(Rcpp::NumericVector y, Rcpp::IntegerVector treatStatus,
                    Rcpp::NumericVector timeRecorded, int iter) {
    TransTest obj(y, treatStatus, timeRecorded);

    obj.runMCMC(iter);

    return Rcpp::List::create(Rcpp::Named("beta0")=obj._beta0,
                              Rcpp::Named("beta1")=obj._beta1,
                              Rcpp::Named("tau")=obj._tau);
}
