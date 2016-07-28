#ifndef _ABC_H_
#define _ABC_H_

#include <vector>
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
        Rcpp::NumericVector _T; // time of recording
        Rcpp::NumericMatrix _D; // data

    public:
        // constructor
        Abc(std::vector<std::vector<double> > data, 
            Rcpp::NumericVector timeRecorded,
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
