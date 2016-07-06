#ifndef _GIBBSSAMPLERTRANSFORMATION_H_
#define _GIBBSSAMPLERTRANSFORMATION_H_

#include "GibbsSamplerMap.h"
#include <vector>
#include <cmath>
#include <Rcpp.h>

class GibbsSamplerTransformation : public GibbsSamplerMap {
  protected:
    // indicates which pattern is growth (will be removed most likely)
    int _whichPattern;
    // vector of case status (assume 0 is non-case) i.e. {0, 0, 0, 1, 1, 1}
    Rcpp::IntegerVector _treatStatus;
    // vector of time recordings
    Rcpp::NumericVector _timeRecorded;

    // abc tolerance
    double _tol;

    // priors for Gibbs Sampling of regression coefficients
    double _mu0;    // beta prior mean
    double _tau0; // beta prior precision
    double _a;      // variance prior shape
    double _b;      // variance prior rate

    // current estimates
    Rcpp::NumericMatrix _beta0;
    Rcpp::NumericMatrix _beta1;
    Rcpp::NumericVector _tau;
    Rcpp::NumericVector _theta;

    // difference
    Rcpp::NumericVector _tolerance;
    Rcpp::NumericVector _proposals;
    Rcpp::NumericVector _accept_prob;

    // logistic growth normalization
    double _normalization;

    // mcmc parameters
    double _delta = 5.0;
    double _prior_mean = 0.0;
    double _prior_sd = 10.0;

  public:
    GibbsSamplerTransformation(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                               double alphaA, double alphaP, double nMaxA, double nMaxP,
                               unsigned long nIterA, unsigned long nIterP,
                               double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                               unsigned long long atomicSize,
                               char label_A, char label_P, char label_D, char label_S,
                               vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                               const string &simulation_id,
                               vector <vector <double> >  &parameters, char the_fixed_matrix, int whichPattern,
                               std::vector<int> treatStatus, std::vector<double> timeRecorded, double tolerance);

    ~GibbsSamplerTransformation() {};

    Rcpp::NumericMatrix beta0();
    Rcpp::NumericMatrix beta1();
    Rcpp::NumericVector theta();
    Rcpp::NumericVector tolerance();
    Rcpp::NumericVector proposals();
    Rcpp::NumericVector accept_prob();
    
    void abc_mcmc(int burn=0, int iter=0, int thin=1, double tolerance=5.0);

};

#endif
