#include "Abc.h"
Abc::Abc(std::vector<std::vector<double> >& data, 
         std::vector<double> timeRecorded,
         Rcpp::NumericVector theta_init,
         std::string prior,
         std::string proposal,
         bool epsilon_mcmc,
         double delta,
         double epsilon,
         double epsilon_rate,
         double prior_mean,
         double prior_sd) :
    _theta(theta_init),
    _pattern(data[0].size()),
    _D(data.size(), timeRecorded.size()) {
    // convert data to Rcpp::NumericMatrix form
    for (unsigned int i = 0; i < data.size(); ++i) {
        for (unsigned int j = 0; j < timeRecorded.size(); ++j) {
            _D(i, j) = data[i][j];
        }
    }

    // initialize _theta to a reasonable value
    // should be parameterized later
    _T=timeRecorded,
    _prior_choice = prior;
    _proposal_choice = proposal;
    _epsilon_mcmc = epsilon_mcmc;
    _delta=delta;
    _epsilon=epsilon;
    _epsilon_rate=epsilon_rate;
    _prior_mean=prior_mean;
    _prior_sd=prior_sd;

    // initialize weights to identity
    old_weight = 1.0;
    new_weight = 1.0;

    // initialize accepted to false
    accepted = false;
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
        return Rcpp::rexp(1, _epsilon_rate)[0];
    } 

    return _epsilon;
}

Rcpp::NumericVector Abc::_epsilon_prior(double param) {
    if (_epsilon_mcmc) {
        Rcpp::NumericVector param_rcpp(1, param);
        return Rcpp::log(Rcpp::dexp(param_rcpp, _epsilon_rate));
    } 

    return Rcpp::wrap(0.0);
}

double Abc::_epsilon_propose() {
    if (_epsilon_mcmc) {
        return Rcpp::rexp(1, _epsilon)[0];
    } 

    return _epsilon;
}

Rcpp::NumericVector Abc::_epsilon_propose(double param1, double param2) {
    if (_epsilon_mcmc) {
        Rcpp::NumericVector param1_rcpp(1, param1);
        return Rcpp::log(Rcpp::dexp(param1_rcpp, param2));
    } 

    return Rcpp::wrap(0.0);
}

void Abc::propose(Rcpp::NumericMatrix A, Rcpp::NumericMatrix P) {

    Rcpp::Rcout << "start empty propose\n";
    // don't update at all
    return;

    // simulate theta' ~ K(theta|theta^{(t-1)})
    Rcpp::NumericVector theta_prime = _proposal();

    // simulate epsilon' ~ K(epsilon|epsilon^{(t-1)})
    double eps_prime = _epsilon_propose();
    eps_prime = std::max(eps_prime, 0.25);

    // simulate x ~ p(x | theta')
    //Rcpp::NumericVector growth = 1. / (1 + Rcpp::exp(-theta_prime * _T));
    Rcpp::NumericMatrix P_prime = P;
    Rcpp::NumericMatrix A_prime = A;

    // allocate curve (assume same number of recordings for all conditions
    int row_num = P.rows();
    P_prime(row_num - 1, Rcpp::_) = curve(theta_prime);

    double Pnorm = Rcpp::sum(P_prime.row(row_num));

    P_prime(row_num - 1, Rcpp::_) = P_prime(row_num - 1, Rcpp::_) / Pnorm;
    A_prime(Rcpp::_, row_num - 1) = A_prime(Rcpp::_, row_num - 1) * Pnorm;

    arma::mat D_prime = Rcpp::as<arma::mat>(A_prime) * Rcpp::as<arma::mat>(P_prime);

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
                 _epsilon_prior(eps_prime) -
                 _epsilon_prior(_epsilon) +
                 _epsilon_propose(_epsilon, eps_prime) -
                 _epsilon_propose(eps_prime, _epsilon) +
                 _proposal(_theta, theta_prime) -
                 _proposal(theta_prime, _theta);
        accept = Rcpp::exp(accept);

        if (u[0] < accept[0]) {
            _theta = theta_prime;
            for (unsigned int i = 0; i < P_prime.cols(); ++i) {
                _pattern[i] = P_prime(2, i);
            }
        } else {
            // c. otherwise
            _theta = _theta;
            //std::fill(_pattern.begin(), _pattern.end(), P(2, Rcpp::_));
            for (unsigned int i = 0; i < P_prime.cols(); ++i) {
                _pattern[i] = P(2, i);
            }
        }
    } else {
    // 4. otherwise
        _theta = _theta;
        for (unsigned int i = 0; i < P_prime.cols(); ++i) {
            _pattern[i] = P(2, i);
        }
    }
}

Rcpp::NumericVector Abc::curve(Rcpp::NumericVector theta_star) {
    Rcpp::NumericVector curve_star(theta_star.length() * _T.length());

    double conds = theta_star.length();
    double recs = _T.length();

    for (int i = 0; i < conds; ++i) {
        for (int j = 0; j < recs; ++j) {
            double tmp = 1.0 / (1.0 + std::exp(-theta_star[i] * _T[j]));
            curve_star[(i + 1) * (j + recs) - recs] = tmp;
        }
    }

    return curve_star;
}

std::vector<double> Abc::pattern() {
    Rcpp::NumericVector tmp = P_update.row(2);
    std::vector<double> tmp2(tmp.begin(), tmp.end());
    return tmp2;
    //return P_update.row(2);
}

std::vector<double> Abc::amplitude() {
    Rcpp::NumericVector tmp = A_update(Rcpp::_, 2);
    std::vector<double> tmp2(tmp.begin(), tmp.end());
    return tmp2;
    // return A_update.col(2);
}

Rcpp::NumericVector Abc::theta() {
    return _theta;
}

Rcpp::NumericVector Abc::epsilon() {
    return _epsilon;
}
