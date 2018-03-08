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
         double prior_sd,
         bool fixed_proposal) :
    _theta(theta_init),
    _theta_truth(theta_init),
    theta_prime(theta_init),
    _D(data.size(), timeRecorded.size() * theta_init.size()) {
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
    _fixed_proposal = fixed_proposal;

    // initialize weights to reflect initial thetas
    old_weight = Rcpp::sum(curve(theta_init));
    new_weight = Rcpp::sum(curve(theta_init));

    _theta_truth = clone(theta_init);
    _theta = clone(theta_init);
    theta_prime = clone(theta_init);

    // initialize accepted to false
    accepted = false;
}

Rcpp::NumericVector Abc::_prior(Rcpp::NumericVector param) {
    Rcpp::NumericVector tmp(1);

    if (_prior_choice == "normal") {
        tmp[0] = Rcpp::sum(Rcpp::dnorm(param, _prior_mean, _prior_sd, true));
    } else if (_prior_choice == "gamma") {
        double a = ((1.0 - _prior_mean) / pow(_prior_sd, 2.0) - 
                    1.0 / _prior_mean) * pow(_prior_mean, 2.0);
        double b = a * (1.0 / _prior_mean - 1.0);
        tmp[0] = Rcpp::sum(Rcpp::dgamma(param, a, b, true));
    } else {
        throw std::logic_error("Invalid prior choice");
    }

    return tmp;
}

Rcpp::NumericVector Abc::_proposal() {
    Rcpp::NumericVector tmp(_theta.length());
    if (_proposal_choice == "normal") {
        for (int i = 0; i < _theta.length(); ++i) {
            tmp[i] = Rcpp::rnorm(1, _theta[i], _delta)[0];
        }
    } else if (_proposal_choice == "gamma") {
        for (int i = 0; i < _theta.length(); ++i) {
            double a = ((1.0 - _theta[i]) / pow(_delta, 2.0) - 
                        1.0 / _theta[i]) * pow(_theta[i], 2.0);
            double b = a * (1.0 / _theta[0] - 1.0);

            tmp[i] = Rcpp::rgamma(1, a, b)[0];
        }
    } else {
        throw std::logic_error("Invalid proposal choice");
    }

    return tmp;
}

Rcpp::NumericVector Abc::_proposal(Rcpp::NumericVector param1,
                                   Rcpp::NumericVector param2) {
    Rcpp::NumericVector tmp(1);
    tmp[0] = 0.0;

    if (_proposal_choice == "normal") {
        // symmetric proposal distribution, so we can ignore it
        // just return additive identity for log scale
    } else if (_proposal_choice == "gamma") {
        for (int i = 0; i < _theta.length(); ++i) {
            double a = ((1.0 - param2[i]) / pow(_delta, 2.0) - 
                        1.0 / param2[i]) * pow(param2[i], 2.0);
            double b = a * (1.0 / param2[i] - 1.0);

            Rcpp::NumericVector arg(1);
            arg[0] = param1[i];

            tmp[0] += Rcpp::dgamma(arg, a, b, true)[0];
        }
    } else {
        throw std::logic_error("Invalid proposal choice");
    }

    return tmp;
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

    accepted = false;

    // simulate theta' ~ K(theta|theta^{(t-1)})
    for (int i = 0; i < _theta_truth.length(); ++i) {
        if (_fixed_proposal) {
            theta_prime[i] = Rcpp::rnorm(1, _theta_truth[i], _delta)[0];
        } else {
            theta_prime[i] = Rcpp::rnorm(1, _theta[i], _delta)[0];
        }
        //theta_prime[i] = 2.91;
    }

    // simulate epsilon' ~ K(epsilon|epsilon^{(t-1)})
    double eps_prime = _epsilon_propose();
    eps_prime = std::max(eps_prime, 0.01);

    // simulate x ~ p(x | theta')
    int row_num = P.rows();
    Rcpp::NumericMatrix P_prime = Rcpp::clone(P);
    Rcpp::NumericMatrix A_prime = Rcpp::clone(A);

    // allocate curve (assume same number of recordings for all conditions
    P_prime(row_num - 1, Rcpp::_) = curve(theta_prime);

    double Pnorm = Rcpp::sum(curve(theta_prime));
    double Pnorm_old;
    if (_fixed_proposal) {
        Pnorm_old = Rcpp::sum(curve(_theta_truth));
    } else {
        Pnorm_old = Rcpp::sum(curve(_theta));
    }

    P_prime(row_num - 1, Rcpp::_) = P_prime(row_num - 1, Rcpp::_) / Pnorm;
    if (weightA) {
        A_prime(Rcpp::_, row_num - 1) = A_prime(Rcpp::_, row_num - 1) * Pnorm / Pnorm_old;
    }

    //Rcpp::Rcout << "Current: " << _theta_truth << " Proposed " << theta_prime << " Pnorm " << Pnorm << " Pnorm_old " << Pnorm_old << "\n";

    arma::mat D_prime = Rcpp::as<arma::mat>(A_prime) * Rcpp::as<arma::mat>(P_prime);
    arma::mat D_orig = Rcpp::as<arma::mat>(A) * Rcpp::as<arma::mat>(P);

    double D_diff = arma::accu(D_prime - D_orig);

    // calculate rho(S(x), S(y))
    arma::mat diff = Rcpp::as<arma::mat>(_D) - D_prime;
    rho = norm(diff, 2);

    // should be calculating eps_prime as some function of
    // current l2norm (i.e. norm(_D - A * P, 2))
    arma::mat orig = Rcpp::as<arma::mat>(_D) - D_orig; 
    rho_thresh = norm(orig, 2) + eps_prime;

    //Rcpp::Rcout << "D_diff " << D_diff << " rho " << rho << " rho_thresh " << rho_thresh << "\n";

    // calculate acceptance probability
    //if (rho < eps_prime) {
    if (rho < rho_thresh) {
        // a. u ~ U(0, 1)
        Rcpp::NumericVector u = Rcpp::runif(1, 0, 1);
        Rcpp::NumericVector accept;
        // TODO: modify accept to vector of length 1

        // b. if u leq pi(theta')/pi*theta^{(t-1)} times 
        //             K(theta^{(t-1)}|theta')/K(theta'|theta^{(t-1)})
        //             theta^{(t)} = theta'
        accept = _prior(theta_prime) -
                 _prior(_theta) +
                 //_epsilon_prior(eps_prime) -
                 //_epsilon_prior(_epsilon) +
                 //_epsilon_propose(_epsilon, eps_prime) -
                 //_epsilon_propose(eps_prime, _epsilon) +
                 _proposal(_theta, theta_prime) -
                 _proposal(theta_prime, _theta);
        accept = Rcpp::exp(accept);

        //Rcpp::Rcout << "acceptance probability: " << accept << "\n";

        if (u[0] < accept[0]) {
            accepted = true;
            _theta = theta_prime;

            //Rcpp::Rcout << "old " << old_weight <<
                           //" new " << new_weight <<
                           //" theta' " << theta_prime <<
                           //" Pnorm " << Pnorm <<
                           //" sum " << Rcpp::sum(curve(theta_prime)) << "\n";

            old_weight = new_weight;
            new_weight = Pnorm;

            //Rcpp::Rcout << "old " << old_weight <<
                           //" new " << new_weight <<
                           //" theta' " << theta_prime <<
                           //" Pnorm " << Pnorm <<
                           //" sum " << Rcpp::sum(curve(theta_prime)) << "\n";
        } else {
            // c. otherwise
            _theta = _theta;
        }
    } else {
    // 4. otherwise
        _theta = _theta;
    }
}

Rcpp::NumericVector Abc::curve(Rcpp::NumericVector theta_star) {
    Rcpp::NumericVector curve_star(theta_star.length() * _T.length());

    double conds = theta_star.length();
    double recs = _T.length();

    for (int i = 0; i < conds; ++i) {
        for (int j = 0; j < recs; ++j) {
            double tmp = 1.0 / (1.0 + std::exp(-theta_star[i] * (_T[j] - 0.5)));
            int ind = i * recs + j;
            curve_star[ind] = tmp;
        }
    }

    return curve_star;
}

Rcpp::NumericVector Abc::theta() {
    return _theta;
}

double Abc::epsilon() {
    return _epsilon;
}
