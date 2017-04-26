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
                                                       std::vector<int> treatStatus, std::vector<double> timeRecorded, 
                                                       std::string prior, std::string proposal,
                                                       bool epsilon_mcmc, double delta, 
                                                       double epsilon, double epsilon_prior, 
                                                       double prior_mean, double prior_sd):
    GibbsSamplerMap(nEquil, nSample, nFactor, alphaA, alphaP, nMaxA, nMaxP, nIterA, nIterP,
                    max_gibbsmass_paraA, max_gibbsmass_paraP, atomicSize, label_A, label_P, label_D, label_S,
                    DVector, SVector, simulation_id, parameters, the_fixed_matrix),
    _theta(nSample + _nEquil),
    _epsilon(nSample + _nEquil),
    _growth(DVector, timeRecorded, prior, proposal, epsilon_mcmc, delta,
            epsilon, epsilon_prior, prior_mean, prior_sd) {

    // assignments for growth data
    _whichPattern = whichPattern;
    _treatStatus = treatStatus;
    _timeRecorded = timeRecorded;

    // store normalization for fixed pattern
    _normalization = 0.0;
    for (int j = 0; j < parameters[0].size(); j++) {
        _normalization += parameters[0][j];
    }
}

Rcpp::NumericVector GibbsSamplerTransformation::theta() {
    return _theta;
}

Rcpp::NumericVector GibbsSamplerTransformation::epsilon() {
    return _epsilon;
}

void GibbsSamplerTransformation::abc_mcmc(int burn, int iter, int thin, double tolerance) {
    Rcpp::Rcout << "allocate A and P\n";
    // get the A, P, D matrices
    Rcpp::NumericMatrix A_curr(_AMatrix.get_nRow(), _AMatrix.get_nCol());
    Rcpp::NumericMatrix P_curr(_PMatrix.get_nRow(), _PMatrix.get_nCol());

    //A
    for (int i = 0; i < A_curr.nrow(); ++i) {
        std::vector<double> curr_row = _AMatrix.get_Row(i);

        for (int j = 0; j < A_curr.ncol(); ++j) {
            A_curr(i, j) = curr_row[j];
        }
    }

    //P
    for (int i = 0; i < P_curr.nrow(); ++i) {
        std::vector<double> curr_row = _PMatrix.get_Row(i);

        for (int j = 0; j < P_curr.ncol(); ++j) {
            P_curr(i, j) = curr_row[j];
        }
    }

    Rcpp::Rcout << "proposal\n";
    for (int s = 0; s < thin; ++s) {
        _growth.propose(A_curr, P_curr);
    }
    Rcpp::Rcout << "Done proposing\n";

    // save the theta's
    _theta[burn + iter] =_growth.theta()[0];
    
    // save the epsilon's
    _epsilon[burn + iter] =_growth.epsilon()[0];

    // update the P matrix
    // _PMatrix.setRow(_growth.pattern(), 2);

    // update the A matrix
    // _AMatrix.setCol(_growth.pattern(), 2);
}
