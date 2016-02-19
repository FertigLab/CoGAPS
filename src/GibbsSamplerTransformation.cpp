#include "GibbsSamplerTransformation.h"

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
                                                       std::vector<int> treatStatus, std::vector<double> timeRecorded) :
    GibbsSamplerMap(nEquil, nSample, nFactor, alphaA, alphaP, nMaxA, nMaxP, nIterA, nIterP,
                    max_gibbsmass_paraA, max_gibbsmass_paraP, atomicSize, label_A, label_P, label_D, label_S,
                    DVector, SVector, simulation_id, parameters, the_fixed_matrix) {
    // assignments for growth data
    _whichPattern = whichPattern;
    _treatStatus = treatStatus;
    _timeRecorded = timeRecorded;

    // initialize priors in constructor
    _mu0 = 0.0;    // beta prior mean
    _tau0 = 0.001; // beta prior precision
    _a = 1.0;      // variance prior shape
    _b = 1.0;      // variance prior rate

    // sample from prior to initialize regression parameters
    int nTreats = Rcpp::unique(treatStatus).size();
    _beta0 = Rcpp::rnorm(nTreats, _mu0, sqrt(1. / _tau0));
    _beta1 = Rcpp::rnorm(nTreats, _mu0, sqrt(1. / _tau0));
    _tau = Rcpp::rgamma(nTreats, _a, _b);
}

void GibbsSamplerTransformation::update_pattern(Rcpp::NumericVector(*transformation)(Rcpp::NumericVector)) {
    // get current pattern data
    Rcpp::NumericVector y_all = _PMatrix.get_Row(_whichPattern);
    int nTreats = Rcpp::unique(_treatStatus).size();

    // http://www.cs.toronto.edu/~radford/csc2541.S11/week3.pdf

    for (int i = 0; i < nTreats; ++i) {
        // initialize lists of y and x for each regression
        Rcpp::NumericVector y = y_all[_treatStatus == i];
        y = transformation(y / (Rcpp::max(y) + 1e-6));
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
        _beta1[i] = Rcpp::as<double>(Rcpp::rnorm(1, post_mean, sqrt(post_var)));

        // \tau | x, y, _beta0, _beta1
        post_shape = _a + n / 2.;
        post_rate = _b + Rcpp::sum(Rcpp::pow(y - _beta0[i] - _beta1[i] * x, 2.0) / 2.0);
        _tau[i] = Rcpp::as<double>(Rcpp::rgamma(1, post_shape, 1. / post_rate));
    }
}
