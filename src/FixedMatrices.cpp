#include "Abc.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix FixedMatrix(Rcpp::NumericMatrix D,
                                Rcpp::NumericMatrix A,
                                Rcpp::NumericMatrix P,
                                Rcpp::NumericVector T,
                                int nSample,
                                double delta,
                                double epsilon,
                                double prior_mean,
                                double prior_sd,
                                bool fixedproposal) {

    std::vector<std::vector<double> > DVector;
    std::vector<double> Tvec = Rcpp::as< std::vector<double> >(T);
    Rcpp::NumericVector theta_init(1, prior_mean);

    int numC = D.ncol();
    int numR = D.nrow();
    double tempFrameElement;
    DVector.resize(numR);

    for (int i = 0; i < numR; i++) {
        DVector[i].resize(numC);
    }

    for (int i = 0; i < numR; i++) {
        for (int j = 0; j < numC; j++) {
            DVector[i][j] = D(i, j);
        }
    }

    Abc growth(DVector, Tvec, theta_init, "normal", "normal", false, 
               delta, epsilon, 0.3, prior_mean, prior_sd);
}

