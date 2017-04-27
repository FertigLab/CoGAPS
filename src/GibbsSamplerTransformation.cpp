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

map <unsigned long long, double> GibbsSamplerTransformation::getADomain() {
    return _AAtomicdomain.getDomain();
}

map <unsigned long long, double> GibbsSamplerTransformation::getPDomain() {
    return _AAtomicdomain.getDomain();
}

map<unsigned int, unsigned long long> GibbsSamplerTransformation::AlBoundariesByBin() {
    return _AAtomicdomain.lBoundariesByBin();
}

map<unsigned long long, unsigned int> GibbsSamplerTransformation::AlBoundaries() {
    return _AAtomicdomain.lBoundaries();
}

map<unsigned int, unsigned long long> GibbsSamplerTransformation::PlBoundariesByBin() {
    return _PAtomicdomain.lBoundariesByBin();
}

map<unsigned long long, unsigned int> GibbsSamplerTransformation::PlBoundaries() {
    return _PAtomicdomain.lBoundaries();
}

void GibbsSamplerTransformation::getAAtomicColumn() {
    // number of rows and columns
    unsigned int rows = _AMatrix.get_nRow();
    unsigned int cols = _AMatrix.get_nCol();

    // get atomic space for A
    map<unsigned long long, double> atoms = _AAtomicdomain.getDomain();
    map<unsigned int, unsigned long long> bins = _AAtomicdomain.lBoundariesByBin();

    // initialize helper variables for accessing space
    unsigned long long lb, ub;
    map<unsigned long long, double>::iterator it = atoms.begin();

    // iterate over the atomic space, but only look at
    // elements mapping to the first column
    for (int i = (cols - 1); i < (rows * cols); i += cols) {
        // track atoms in cell of last column
        int hits = 0;

        // show row number
        Rcpp::Rcout << (i - cols + 1) / 3 << " ";

        // bounds on bin
        lb = bins[i];
        if (i == (rows * cols - 1)) {
            ub = std::numeric_limits<unsigned long long>::max();
        } else {
            ub = bins[i+1];
        }

        // get the mass of the atoms within the bin
        while (it->first < ub & it != atoms.end()) {
            if (it->first >= lb) {
                hits++; // found an atom
                Rcpp::Rcout << it->second << " "; // print mass
            }

            it++;
        }

        // if not atoms, print a 0
        if (hits == 0) {
            Rcpp::Rcout << "0.0";
        }

        Rcpp::Rcout << "\n";
    }
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
