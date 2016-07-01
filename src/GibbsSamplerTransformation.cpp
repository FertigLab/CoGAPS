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
    _theta(nSample + _nEquil),
    _tolerance(nSample + _nEquil),
    _proposals(nSample),
    _accept_prob(nSample + _nEquil) {
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

    _theta[0] = Rcpp::as<double>(Rcpp::rnorm(1, _prior_mean, _prior_sd));
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

Rcpp::NumericVector GibbsSamplerTransformation::tolerance() {
    return _tolerance;
}

Rcpp::NumericVector GibbsSamplerTransformation::proposals() {
    return _proposals;
}

Rcpp::NumericVector GibbsSamplerTransformation::accept_prob() {
    return _accept_prob;
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
    _proposals[iter] = theta1;

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

    _tolerance[iter] = d1;

    if (d1 < _tol) {
        _theta[iter] = theta1;
        // accept
        _PMatrix.setRow(Rcpp::as<std::vector<double> >(logit_patt_1), 2);
    } 
    // else reject
}

void GibbsSamplerTransformation::update_pattern_abc_mcmc(int burn, int iter) {
    // get the A, P, D matrices
    Rcpp::NumericMatrix A_curr(_AMatrix.get_nRow(), _AMatrix.get_nCol());
    Rcpp::NumericMatrix P_curr(_PMatrix.get_nRow(), _PMatrix.get_nCol());
    Rcpp::NumericMatrix D(_DMatrix.get_nRow(), _DMatrix.get_nCol());

    //A
    for (int i = 0; i < A_curr.nrow(); ++i) {
        std::vector<double> curr_row = _AMatrix.get_Row(i);

        for (int j = 0; j < A_curr.ncol(); ++j) {
            A_curr(i, j) = curr_row[j];
        }
    }

    //P
    for (int i = 0; i < P_curr.nrow(); ++i) {
        std::vector<double> curr_row = _PMatrix.get_Row(i);

        for (int j = 0; j < P_curr.ncol(); ++j) {
            P_curr(i, j) = curr_row[j];
        }
    }

    // D
    for (int i = 0; i < D.nrow(); ++i) {
        std::vector<double> curr_row = _DMatrix.get_Row(i);

        for (int j = 0; j < D.ncol(); ++j) {
            D(i, j) = curr_row[j];
        }
    }

    // last theta
    Rcpp::NumericVector theta = _theta[burn + iter - 1];

    // 1. simulate theta' ~ K(theta|theta^{(t-1)})
    Rcpp::NumericVector theta_prime = Rcpp::rnorm(1, theta[0], _delta);
    _proposals[burn + iter] = theta_prime[0];

    // 2. simulate x ~ p(x | theta')
    Rcpp::NumericVector T = _timeRecorded[_treatStatus == 0];
    Rcpp::NumericVector growth = 1. / (1 + Rcpp::exp(-theta_prime * T));
    Rcpp::NumericMatrix P_prime = P_curr;

    for (int i = 0; i < 10; ++i) {
        P_prime(2, i) = growth[i];
    }

    // 3. If rho(S(x), S(y)) < epsilon
    arma::mat D_prime = Rcpp::as<arma::mat>(A_curr) * Rcpp::as<arma::mat>(P_prime);
    arma::mat diff = Rcpp::as<arma::mat>(D) - D_prime;
    double rho = norm(diff, 2);

    _tolerance[burn + iter] = rho;

    if (rho < _tol) {
        // a. u ~ U(0, 1)
        Rcpp::NumericVector u = Rcpp::runif(1, 0, 1);

        // b. if u leq pi(theta')/pi*theta^{(t-1)} times 
        //             K(theta^{(t-1)}|theta')/K(theta'|theta^{(t-1)})
        //             theta^{(t)} = theta'
        Rcpp::NumericVector accept = Rcpp::dnorm(theta_prime, _prior_mean, _prior_sd, true) -
                                     Rcpp::dnorm(theta, _prior_mean, _prior_sd, true) +
                                     Rcpp::dnorm(theta, theta_prime[0], _delta) -
                                     Rcpp::dnorm(theta_prime, theta[0], _delta);
        accept = Rcpp::exp(accept);
        _accept_prob[burn + iter] = accept[0];

        if (u[0] < accept[0]) {
            _theta[burn + iter] = theta_prime[0];
        } else {
            // c. otherwise
            _theta[burn + iter] = theta[0];
        }
    } else {
    // 4. otherwise
        _theta[burn + iter] = theta[0];
    }
}
