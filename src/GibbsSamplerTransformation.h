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

    // priors for Gibbs Sampling of regression coefficients
    double _mu0;    // beta prior mean
    double _tau0; // beta prior precision
    double _a;      // variance prior shape
    double _b;      // variance prior rate

    // current estimates
    Rcpp::NumericMatrix _beta0;
    Rcpp::NumericMatrix _beta1;
    Rcpp::NumericVector _tau;

    // logistic growth normalization
    double _normalization;

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
                               std::vector<int> treatStatus, std::vector<double> timeRecorded);

    ~GibbsSamplerTransformation() {};

    Rcpp::NumericMatrix beta0();
    Rcpp::NumericMatrix beta1();
    
    // expect that some transformation of the pattern will be linear so
    // allow for passing an a priori function that transforms the pattern
    // to a linear pattern
    // for example, if we expect the pattern to be logistic, we can pass a logistic
    // function to transform the data i.e.
    // GibbSampPatt.update_pattern(GibbsSampPatt::&logit);
    void update_pattern(Rcpp::NumericVector(*transformation)(Rcpp::NumericVector), int iter=0);

    void update_pattern_test(Rcpp::NumericVector(*transformation)(Rcpp::NumericVector), int iter=0);
    void update_pattern_abc(Rcpp::NumericVector(*transformation)(Rcpp::NumericVector), int iter=0);

};

#endif
