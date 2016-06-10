#include "GibbsSamplerTransformation.h"
#include <Rcpp.h>
#include <RcppArmadillo.h>

// ******************** CONSTRUCTOR ********************************************
GibbsSamplerTransformation::GibbsSamplerTransformation(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                                                       double alphaA, double alphaP, double nMaxA, double nMaxP,
                                                       unsigned long nIterA, unsigned long nIterP,
                                                       double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                                                       unsigned long long atomicSize,
                                                       char label_A, char label_P, char label_D, char label_S,
                                                       vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                                                       const string &simulation_id,
                                                       vector <vector <double> >  &parameters, char the_fixed_matrix, int whichPattern,
                                                       std::vector<int> treatStatus, std::vector<double> timeRecorded, double tolerance) :
    GibbsSamplerMap(nEquil, nSample, nFactor, alphaA, alphaP, nMaxA, nMaxP, nIterA, nIterP,
                    max_gibbsmass_paraA, max_gibbsmass_paraP, atomicSize, label_A, label_P, label_D, label_S,
                    DVector, SVector, simulation_id, parameters, the_fixed_matrix),
    _beta0(nSample, 2),
    _beta1(nSample, 2),
    _theta(nSample) {
    // assignments for growth data
    _whichPattern = whichPattern;
    _treatStatus = treatStatus;
    _timeRecorded = timeRecorded;

    // store normalization for fixed pattern
    _normalization = 0.0;
    for (int j = 0; j < parameters[0].size(); j++) {
        _normalization += parameters[0][j];
    }

    // abc tolerance
    _tol = tolerance;

    // initialize priors in constructor
    _mu0 = 0.0;    // beta prior mean
    _tau0 = 0.001; // beta prior precision
    _a = 0.001;      // variance prior shape
    _b = 0.001;      // variance prior rate

    // sample from prior to initialize regression parameters
    int nTreats = Rcpp::unique(_treatStatus).size();
    
    _beta0(0, Rcpp::_) = Rcpp::rnorm(nTreats, _mu0, sqrt(1. / _tau0));
    _beta1(0, Rcpp::_) = Rcpp::rnorm(nTreats, _mu0, sqrt(1. / _tau0));
    _tau = Rcpp::rgamma(nTreats, _a, _b);
}

Rcpp::NumericMatrix GibbsSamplerTransformation::beta0() {
    return _beta0;
}

Rcpp::NumericMatrix GibbsSamplerTransformation::beta1() {
    return _beta1;
}

Rcpp::NumericVector GibbsSamplerTransformation::theta() {
    return _theta;
}

void GibbsSamplerTransformation::update_pattern(Rcpp::NumericVector(*transformation)(Rcpp::NumericVector),
                                                int iter) {
    // get current pattern data
    // maps stores fixed pattern starting from last row. Let's assume that theres only one fixed pattern
    std::vector<double> y_all_temp = _PMatrix.get_Row(_nFactor - 1);
    Rcpp::NumericVector y_all;
    y_all = y_all_temp;
    int nTreats = Rcpp::unique(_treatStatus).size();

    // since iter defaults to zero, create a past iteration variable
    int past_iter = 0;
    if (iter > past_iter) {
        past_iter = iter - 1;
    }

    // http://www.cs.toronto.edu/~radford/csc2541.S11/week3.pdf

    for (int i = 0; i < nTreats; ++i) {
        // initialize lists of y and x for each regression
        Rcpp::NumericVector y = y_all[_treatStatus == i];
        // "de-normalize" data
        y = y *  _normalization;
        //y = transformation(y / (Rcpp::max(y) + 1e-16));
        y = transformation(y);
        Rcpp::NumericVector x = _timeRecorded[_treatStatus == i];
        Rcpp::NumericVector x_sq = Rcpp::pow(x, 2.0);

        // initialize variables for full conditionals
        double post_mean, post_var;     // normal distributions
        double post_shape, post_rate;   // gamma distributions
        int n = y.size();

        // \beta_0 | x, y, \beta_1, \tau
        post_var = 1. / (_tau0 + n * _tau[i]);
        post_mean = (_tau0 * _mu0 + _tau[i] * Rcpp::sum(y - _beta1(past_iter, i) * x)) * post_var;
        _beta0(iter, i) = Rcpp::as<double>(Rcpp::rnorm(1, post_mean, 
                                                 sqrt(post_var)));
        // \beta_1 | x, y, \beta_0, \tau
        post_var = 1. / (_tau0 + _tau[i] * Rcpp::sum(x_sq));
        post_mean = (_tau0 * _mu0 + _tau[i] * Rcpp::sum(x * (y - _beta1(past_iter, i)))) * post_var;
        _beta1(iter, i) = Rcpp::as<double>(Rcpp::rnorm(1, post_mean, sqrt(post_var)));

        // \tau | x, y, _beta0, _beta1
        post_shape = _a + n / 2.;
        post_rate = _b + Rcpp::sum(Rcpp::pow(y - _beta0(past_iter, i) - _beta1(past_iter, i) * x, 2.0) / 2.0);
        _tau[i] = Rcpp::as<double>(Rcpp::rgamma(1, post_shape, 1. / post_rate));
    }
}

void GibbsSamplerTransformation::update_pattern_abc(Rcpp::NumericVector(*transformation)(Rcpp::NumericVector),
                                                    int iter) {
    // since iter defaults to zero, create a past iteration variable
    int past_iter = 0;
    if (iter > past_iter) {
        past_iter = iter - 1;
    }

    // propose new parameters
    double theta1;
    theta1 = std::abs(Rcpp::as<double>(Rcpp::rnorm(1, 0, 10)));

    // now build logistic growth
    Rcpp::NumericVector x = _timeRecorded[_treatStatus == 0];
    Rcpp::NumericVector patt1 = 1. / (1 + Rcpp::exp(-theta1 * x));

    // next let's get the current A and P matrices
    Rcpp::NumericMatrix A_curr(_AMatrix.get_nRow(), _AMatrix.get_nCol());
    Rcpp::NumericMatrix P_curr(_PMatrix.get_nRow(), _PMatrix.get_nCol());

    for (int i = 0; i < A_curr.nrow(); ++i) {
        std::vector<double> curr_row = _AMatrix.get_Row(i);

        for (int j = 0; j < A_curr.ncol(); ++j) {
            A_curr(i, j) = curr_row[j];
        }
    }

    for (int i = 0; i < P_curr.nrow(); ++i) {
        std::vector<double> curr_row = _PMatrix.get_Row(i);

        for (int j = 0; j < P_curr.ncol(); ++j) {
            P_curr(i, j) = curr_row[j];
        }
    }

    // construct data matrix
    Rcpp::NumericMatrix D(_DMatrix.get_nRow(), _DMatrix.get_nCol());

    for (int i = 0; i < D.nrow(); ++i) {
        std::vector<double> curr_row = _DMatrix.get_Row(i);

        for (int j = 0; j < D.ncol(); ++j) {
            D(i, j) = curr_row[j];
        }
    }

    // replace last row of P_curr, one pattern at a time
    Rcpp::NumericMatrix P1 = P_curr;

    // SEGFAULT IS HERE
    // get elements of past logit pattern
    Rcpp::NumericVector logit_patt_1 = Rcpp::wrap(_PMatrix.get_Row(P_curr.nrow() - 1));

    // plug in current simulated values
    for (int i = 0; i < 10; ++i) {
        logit_patt_1[i] = patt1[i];
    }

    // now put the matrices back together
    P1(2, Rcpp::_) = logit_patt_1;

    // now perform matrix multiplication
    arma::mat D_prime1 = Rcpp::as<arma::mat>(A_curr) * Rcpp::as<arma::mat>(P1);

    // calculate summary statistics/distance

    arma::mat D_diff1 = Rcpp::as<arma::mat>(D) - D_prime1;

    double d1;
    d1 = norm(D_diff1, 2);

    if (d1 < _tol) {
        _theta[iter] = theta1;
        // accept
        _PMatrix.setRow(Rcpp::as<std::vector<double> >(logit_patt_1), 2);
    } 
    // else reject
}
