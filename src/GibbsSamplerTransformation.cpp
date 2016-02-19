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

void GibbsSamplerTransformation::update_pattern(std::vector<double>(*transformation)(std::vector<double>)) {
    // get current pattern data
    arma::vec y_all = _PMatrix.get_Row(_whichPattern);

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
        _beta0[i] = randgen('N', post_mean, sqrt(post_var));

        // \beta_1 | x, y, \beta_0, \tau
        post_var = 1. / (_tau0 + _tau[i] * arma::sum(x_sq));
        post_mean = (_tau0 * _mu0 + _tau[i] * arma::sum(x * (y - _beta1[i]))) * post_var;
        _beta1[i] = randgen('N', post_mean, sqrt(post_var));

        // \tau | x, y, _beta0, _beta1
        post_shape = _a + n / 2.;
        post_rate = _b + arma::sum(arma::pow(y - _beta0[i] - _beta1[i] * x, 2.0) / 2.0);
        _tau[i] = arma::conv_to<double>::from(arma::randg(1, arma::distr_param(post_shape, post_rate)));
    }
}
