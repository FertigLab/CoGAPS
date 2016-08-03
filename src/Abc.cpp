#include "Abc.h"
Abc::Abc(std::vector<std::vector<double> >& data, 
         std::vector<double> timeRecorded,
         std::string prior,
         std::string proposal,
         bool epsilon_mcmc,
         double delta,
         double epsilon,
         double epsilon_rate,
         double prior_mean,
         double prior_sd) :
    _theta(1),
    _D(data.size(), timeRecorded.size()) {
    // convert data to Rcpp::NumericMatrix form
    for (unsigned int i = 0; i < data.size(); ++i) {
        for (unsigned int j = 0; j < timeRecorded.size(); ++j) {
            _D(i, j) = data[i][j];
        }
    }

    _T=timeRecorded,
    _prior_choice = prior;
    _proposal_choice = proposal;
    _epsilon_mcmc = epsilon_mcmc;
    _delta=delta;
    _epsilon=epsilon;
    _epsilon_rate=epsilon_rate;
    _prior_mean=prior_mean;
    _prior_sd=prior_sd;
}

Rcpp::NumericVector Abc::_prior(Rcpp::NumericVector param) {
    if (_prior_choice == "normal") {
        return Rcpp::dnorm(param, _prior_mean, _prior_sd, true);
    } else if (_prior_choice == "gamma") {
        double a = ((1.0 - _prior_mean) / pow(_prior_sd, 2.0) - 
                    1.0 / _prior_mean) * pow(_prior_mean, 2.0);
        double b = a * (1.0 / _prior_mean - 1.0);
        return Rcpp::dgamma(param, a, b, true);
    } else {
        throw std::logic_error("Invalid prior choice");
        return 1;
    }
}

Rcpp::NumericVector Abc::_proposal() {
    if (_proposal_choice == "normal") {
        return Rcpp::rnorm(1, _theta[0], _delta);
    } else if (_proposal_choice == "gamma") {
        double a = ((1.0 - _theta[0]) / pow(_delta, 2.0) - 
                    1.0 / _theta[0]) * pow(_theta[0], 2.0);
        double b = a * (1.0 / _theta[0] - 1.0);
        return Rcpp::rgamma(1, a, b);
    } else {
        throw std::logic_error("Invalid proposal choice");
        return 1;
    }
}

Rcpp::NumericVector Abc::_proposal(Rcpp::NumericVector param1,
                              Rcpp::NumericVector param2) {
    if (_proposal_choice == "normal") {
        // symmetric proposal distribution, so we can ignore it
        // just return additive identity for log scale
        return Rcpp::wrap(0);
    } else if (_proposal_choice == "gamma") {
        double a = ((1.0 - param2[0]) / pow(_delta, 2.0) - 
                    1.0 / param2[0]) * pow(param2[0], 2.0);
        double b = a * (1.0 / param2[0] - 1.0);
        return Rcpp::dgamma(param1, a, b, true);
    } else {
        throw std::logic_error("Invalid proposal choice");
        return 1;
    }
}

double Abc::_epsilon_prior() {
    if (_epsilon_mcmc) {
        return Rcpp::rexp(1, _epsilon_rate);
    } 

    return _epsilon;
}

Rcpp::NumericVector Abc::_epsilon_prior(double param) {
    if (_epsilon_mcmc) {
        return Rcpp::dexp(param, _epsilon_rate, log=true);
    } 

    return Rcpp::wrap(0.0);
}

double Abc::_epsilon_propose() {
    if (_epsilon_mcmc) {
        return Rcpp::rexp(1, 1 / _epsilon)[0];
    } 

    return _epsilon;
}

Rcpp::NumericVector Abc::_epsilon_propose(double param1, double param2) {
    if (_epsilon_mcmc) {
        Rcpp::dexp(Rcpp::wrap(param1), 1.0 / param2, log=true)
    } 

    return Rcpp::wrap(0.0);
}

void Abc::propose(Rcpp::NumericMatrix A, Rcpp::NumericMatrix P) {
    // simulate theta' ~ K(theta|theta^{(t-1)})
    Rcpp::NumericVector theta_prime = _proposal();

    // simulate epsilon' ~ K(epsilon|epsilon^{(t-1)})
    double eps_prime = _epsilon_propose();
    eps_prime = std::max(eps_prime, 0.25);

    // simulate x ~ p(x | theta')
    Rcpp::NumericVector growth = 1. / (1 + Rcpp::exp(-theta_prime * _T));
    Rcpp::NumericMatrix P_prime = P;

    for (int i = 0; i < 10; ++i) {
        P_prime(2, i) = growth[i];
    }

    arma::mat D_prime = Rcpp::as<arma::mat>(A) * Rcpp::as<arma::mat>(P_prime);

    // calculate rho(S(x), S(y))
    arma::mat diff = Rcpp::as<arma::mat>(_D) - D_prime;
    double rho = norm(diff, 2);

    // calculate acceptance probability
    if (rho < eps_prime) {
        // a. u ~ U(0, 1)
        Rcpp::NumericVector u = Rcpp::runif(1, 0, 1);
        Rcpp::NumericVector accept;

        // b. if u leq pi(theta')/pi*theta^{(t-1)} times 
        //             K(theta^{(t-1)}|theta')/K(theta'|theta^{(t-1)})
        //             theta^{(t)} = theta'
        accept = _prior(theta_prime) -
                 _prior(_theta) +
                 _proposal(_theta, theta_prime) -
                 _proposal(theta_prime, _theta);
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

std::vector<double> Abc::pattern() {
    Rcpp::NumericVector growth = 1. / (1 + Rcpp::exp(_theta * _T));

    return Rcpp::as<std::vector<double> >(growth);
}

Rcpp::NumericVector Abc::theta() {
    return _theta;
}

Rcpp::NumericVector Abc::epsilon() {
    return _epsilon;
}
