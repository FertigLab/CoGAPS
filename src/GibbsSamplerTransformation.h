#ifndef _GIBBSSAMPLERTRANSFORMATION_H_
#define _GIBBSSAMPLERTRANSFORMATION_H_

#include "GibbsSamplerMap.h"
#include <vector>
#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

class GibbsSamplerTransformation : public GibbsSamplerMap {
  protected:
    // indicates which pattern is growth (will be removed most likely)
    int _whichPattern;
    // vector of case status (assume 0 is non-case) i.e. {0, 0, 0, 1, 1, 1}
    std::vector<int> _treatStatus;
    // vector of time recordings
    std::vector<double> _timeRecorded;

    // priors for Gibbs Sampling of regression coefficients
    double _mu0 = 0.0;    // beta prior mean
    double _tau0 = 0.0;   // beta prior precision
    double _a = 1.0;      // variance prior shape
    double _b = 1.0;      // variance prior rate

    // current beta estimates
    arma::vec _beta0;
    arma::vec _beta1;

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
                               std::vector<int> caseStatus, std::vector<double> timeRecorded);

    ~GibbsSamplerTransformation() {};
    
    // expect that some transformation of the pattern will be linear so
    // allow for passing an a priori function that transforms the pattern
    // to a linear pattern
    // for example, if we expect the pattern to be logistic, we can pass a logistic
    // function to transform the data i.e.
    // GibbSampPatt.update_pattern(GibbsSampPatt::&logit);
    void update_pattern(std::vector<double>(*transformation)(std::vector<double>));

    // transformations
    // logistic transformation (log(p / 1 - p))
    std::vector<double> logit(std::vector<double> data);
    // identity transformation (pattern is already linear)
    std::vector<double> identity(std::vector<double> data);
};

#endif
