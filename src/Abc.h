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

        // data
        Rcpp::NumericVector _T;
        Rcpp::NumericMatrix _D;

    public:
        // propose new theta
        void propose(Rcpp::NumericMatrix A, Rcpp::NumericMatrix P);

        // get theta
        Rcpp::NumericVector theta();

};

#endif
