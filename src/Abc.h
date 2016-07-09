#ifndef _ABC_H_
#define _ABC_H_

#include <Rcpp.h>

class Abc {
    private:
        // MH parameters
        double _delta;

        // ABC parameters
        double _epsilon;

        // value 
        Rcpp::NumericVector _theta;

        // data
        Rcpp::NumericVector _timeRecorded;
        Rcpp::NumericVector _treatStatus;
        Rcpp::NumericMatrix _D;

    public:
        // propose new theta
        void propose(Rcpp::NumericMatrix A, Rcpp::NumericMatrix P);

        // get theta
        Rcpp::NumericVector theta();

};

#endif
