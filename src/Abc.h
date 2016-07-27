#ifndef _ABC_H_
#define _ABC_H_

#include <Rcpp.h>
#include <RcppArmadillo.h>

class Abc {
    private:
        // MH parameters
        double _delta;

        // ABC parameters
        double _epsilon;

        // value 
        Rcpp::NumericVector _theta;

        // priors
        double _prior_mean;
        double _prior_sd;

        // data
        // Rcpp::NumericVector _T;
        Rcpp::NumericMatrix _D;

    public:
        // constructor
        Abc::Abc(Rcpp::NumericMatrix data, 
                 double delta=10.0,
                 double epsilon=100.0,
                 double prior_mean=0.0,
                 double prior_sd=10.0);

        // propose new theta
        void propose(Rcpp::NumericMatrix A, Rcpp::NumericMatrix P);

        // get parameters
        Rcpp::NumericVector theta();
        Rcpp::NumericVector epsilon();

};

#endif
