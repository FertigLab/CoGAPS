#include "Abc.h"
Abc::Abc(std::vector<std::vector<double> >& data, 
         Rcpp::NumericVector timeRecorded,
         double delta,
         double epsilon,
         double prior_mean,
         double prior_sd) :
    _theta(1),
    _D(data.size(), timeRecorded.size()) {
    // convert data to Rcpp::NumericMatrix form
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < timeRecorded.size(); ++j) {
            _D(i, j) = data[i][j];
        }
    }

    _T=timeRecorded,
    _delta=delta;
    _epsilon=epsilon;
    _prior_mean=prior_mean;
    _prior_sd=prior_sd;
}

void Abc::propose(Rcpp::NumericMatrix A, Rcpp::NumericMatrix P) {
    // 1. simulate theta' ~ K(theta|theta^{(t-1)})
    Rcpp::NumericVector theta_prime = Rcpp::rnorm(1, _theta[0], _delta);

    // 2. simulate x ~ p(x | theta')
    Rcpp::NumericVector growth = 1. / (1 + Rcpp::exp(-theta_prime * _T));
    Rcpp::NumericMatrix P_prime = P;

    for (int i = 0; i < 10; ++i) {
        P_prime(2, i) = growth[i];
    }

    // 3. If rho(S(x), S(y)) < epsilon
    arma::mat D_prime = Rcpp::as<arma::mat>(A) * Rcpp::as<arma::mat>(P_prime);
    arma::mat diff = Rcpp::as<arma::mat>(_D) - D_prime;
    double rho = norm(diff, 2);

    if (rho < _epsilon) {
        // a. u ~ U(0, 1)
        Rcpp::NumericVector u = Rcpp::runif(1, 0, 1);
        Rcpp::NumericVector accept;

        // b. if u leq pi(theta')/pi*theta^{(t-1)} times 
        //             K(theta^{(t-1)}|theta')/K(theta'|theta^{(t-1)})
        //             theta^{(t)} = theta'
        accept = Rcpp::dnorm(theta_prime, _prior_mean, _prior_sd, true) -
                 Rcpp::dnorm(_theta, _prior_mean, _prior_sd, true) +
                 Rcpp::dnorm(_theta, theta_prime[0], _delta) -
                 Rcpp::dnorm(theta_prime, _theta[0], _delta);
        accept = Rcpp::exp(accept);

        if (u[0] < accept[0]) {
            _theta = theta_prime;
        } else {
            // c. otherwise
            _theta = _theta;
        }
    } else {
    // 4. otherwise
        _theta = _theta;
    }
}

Rcpp::NumericVector Abc::theta() {
    return _theta;
}

Rcpp::NumericVector Abc::epsilon() {
    return _epsilon;
}
