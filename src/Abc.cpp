#include "Abc.h"

void Abc::propose(Rcpp::NumericMatrix A, Rcpp::NumericMatrix P) {
    // 1. simulate theta' ~ K(theta|theta^{(t-1)})
    Rcpp::NumericVector theta_prime = Rcpp::rnorm(1, _theta[0], _delta);

    // 2. simulate x ~ p(x | theta')
    Rcpp::NumericVector growth = 1. / (1 + Rcpp::exp(-theta_prime * _T));
    Rcpp::NumericMatrix P_prime = P;

    for (int i = 0; i < 10; ++i) {
        P_prime(2, i) = growth[i];
    }

}
