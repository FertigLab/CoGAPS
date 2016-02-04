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

    // size regression parameters
    _beta0.resize(nFactor);
    _beta1.resize(nFactor);
    _tau.resize(nFactor);

    // sample from prior to initialize regression parameters
    for (int i = 0; i < nFactor; ++i) {
        _beta0[i] = randgen("N", _mu0, sqrt(1. / _tau0));
        _beta1[i] = randgen("N", _mu0, sqrt(1. / _tau0));
        _tau[i] = arma::randg(1, arma::distr_param(_a, _b));
    }
}

std::vector<double> GibbsSamplerTransformation::logit(std::vector<double> data) {
    std::vector<double> transformation(data.size());

    for (int i = 0; i < data.size(); ++i) {
        transformation[i] = log(data[i] / (1 - data[i]));
    }

    return transformation;
}

std::vector<double> GibbsSamplerTransformation::identity(std::vector<double> data) {
    return data;
}

void GibbsSamplerTransformation::update_pattern(std::vector<double>(*transformation)(std::vector<double>)) {
    // get current pattern data
    double** POrig = _PMatrix.get_matrix();
    arma::vec y_all = POrig[_whichPattern, ];

    // split data into lists
    for (int i = 0; i < _nFactor; ++i) {
        _y[i] = y.elem(y_all.find(_treatStatus == i));
        _x[i] = _timeRecorded.elem(_timeRecorded.find(_treatStatus == i));
    }

    // split to do ....

    for (int i = 0; i < pattern_by_case.size(); ++i) {
        // initialize variables for full conditionals
        double post_mean, post_var;     // normal distributions
        double post_shape, post_rate;   // gamma distributions
        int n = _y[i].size();

        // \beta_0 | x, y, \beta_1, \tau
        post_var = 1. / (_tau0 + n * _tau);
        post_mean = (_tau0 * _mu0 + _tau * arma::sum(_y[i] - _beta1 * _x[i])) * post_var;
        _beta0[i] = randgen("N", post_mean, sqrt(post_var));

        // \beta_1 | x, y, \beta_0, \tau
        post_var = 1. / (_tau0 + _tau * arma::sum(pow(_x[i], 2.0)))
        post_mean = (_tau0 * _mu0 + _tau * arma::sum(_x[i] * (_y[i] - _beta1)) * post_var;
        _beta1[i] = randgen("N", post_mean, sqrt(post_var));

        // \tau | x, y, _beta0, _beta1
        post_shape = _a + n / 2.;
        post_rate = _b + arma::sum(pow(_y[i] - _beta0 - _beta1 * _x[i], 2.0) / 2.0);
        _tau[i] = arma::randg(1, arma::distr_param(post_shape, post_rate));
    }
}
