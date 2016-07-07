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

    // mcmc parameters
    _delta = 5.0;
    _prior_mean = 0.0;
    _prior_sd = 10.0;

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

void GibbsSamplerTransformation::abc_mcmc(int burn, int iter, int thin, double tolerance) {
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
    Rcpp::NumericVector theta(1);
    theta[0] = _theta[max(burn + iter - 1, 0)];

    // initialize mcmc params
    Rcpp::NumericVector theta_prime;
    Rcpp::NumericVector accept; 

    for (int s = 0; s < thin; ++s) {
        // 1. simulate theta' ~ K(theta|theta^{(t-1)})
        theta_prime = Rcpp::rnorm(1, theta[0], _delta);

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
            accept = Rcpp::dnorm(theta_prime, _prior_mean, _prior_sd, true) -
                     Rcpp::dnorm(theta, _prior_mean, _prior_sd, true) +
                     Rcpp::dnorm(theta, theta_prime[0], _delta) -
                     Rcpp::dnorm(theta_prime, theta[0], _delta);
            accept = Rcpp::exp(accept);

            if (u[0] < accept[0]) {
                theta = theta_prime;
            } else {
                // c. otherwise
                theta = theta;
            }
        } else {
        // 4. otherwise
            theta = theta;
        }
    }

    _theta[burn + iter] = theta[0];
    _accept_prob[burn + iter] = accept[0];
    _proposals[burn + iter] = theta_prime[0];
}
