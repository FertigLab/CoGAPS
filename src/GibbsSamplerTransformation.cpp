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
                                                       Rcpp::NumericVector theta_init,
                                                       std::string prior, std::string proposal,
                                                       bool epsilon_mcmc, double delta, 
                                                       double epsilon, double epsilon_prior, 
                                                       double prior_mean, double prior_sd,
                                                       bool fixedproposal):
    GibbsSamplerMap(nEquil, nSample, nFactor, alphaA, alphaP, nMaxA, nMaxP, nIterA, nIterP,
                    max_gibbsmass_paraA, max_gibbsmass_paraP, atomicSize, label_A, label_P, label_D, label_S,
                    DVector, SVector, simulation_id, parameters, the_fixed_matrix),
    _theta(nSample + _nEquil, theta_init.length()),
    _thresh(nSample + _nEquil, 2),
    _epsilon(nSample + _nEquil),
    _growth(DVector, timeRecorded, theta_init, prior, proposal, epsilon_mcmc, delta,
            epsilon, epsilon_prior, prior_mean, prior_sd, fixedproposal) {

    // assignments for growth data
    _whichPattern = whichPattern;
    _treatStatus = treatStatus;
    _timeRecorded = timeRecorded;

    // store normalization for fixed pattern
    _normalization = 0.0;
    for (unsigned int j = 0; j < parameters[0].size(); j++) {
        _normalization += parameters[0][j];
    }

    accepted = 0;
    proposed = 0;
}

Rcpp::NumericMatrix GibbsSamplerTransformation::theta() {
    return _theta;
}

Rcpp::NumericMatrix GibbsSamplerTransformation::thresh() {
    return _thresh;
}

Rcpp::NumericVector GibbsSamplerTransformation::epsilon() {
    return _epsilon;
}

map <unsigned long long, double> GibbsSamplerTransformation::getADomain() {
    return _AAtomicdomain.getDomain();
}

map <unsigned long long, double> GibbsSamplerTransformation::getPDomain() {
    return _PAtomicdomain.getDomain();
}

void GibbsSamplerTransformation::setMass(unsigned long long location, double weight) {
    _AAtomicdomain.setMass(location, weight);
}

double GibbsSamplerTransformation::calcWeight() {
    return _growth.new_weight / _growth.old_weight;
}

void GibbsSamplerTransformation::weightAColumn(double weight) {
    if (_growth.accepted == true) {
        unsigned int col = _AMatrix.get_nCol() - 1;
        unsigned int row = _AMatrix.get_nRow();
        std::vector<double> col_val = _AMatrix.get_Col(col);

        for (unsigned int i = 0; i < row; ++i) {
            col_val[i] *= weight;
        }

        _AMatrix.setCol(col_val, col);
    }
}

void GibbsSamplerTransformation::updatePRow() {
    if (_growth.accepted == true) {
        Rcpp::NumericVector tmp = _growth.curve(_growth.theta()) / _growth.new_weight;
        std::vector<double> new_curve = Rcpp::as<std::vector<double> >(tmp);
        _PMatrix.setRow(new_curve, _PMatrix.get_nRow()-1);
    }
}

void GibbsSamplerTransformation::weightAAtomicColumn(double weight) {
    if (_growth.accepted == true) {
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
        // elements mapping to the last column
        for (unsigned int i = (cols - 1); i < (rows * cols); i += cols) {
            // bounds on bin
            lb = bins[i];
            if (i == (rows * cols - 1)) {
                ub = std::numeric_limits<unsigned long long>::max();
            } else {
                ub = bins[i+1];
            }

            // get the mass of the atoms within the bin
            while ((it->first < ub) & (it != atoms.end())) {
                if (it->first >= lb) {
                    // reweight mass
                    _AAtomicdomain.setMass(it->first, weight);
                }

                it++;
            }
        }
    }
}

Rcpp::NumericMatrix GibbsSamplerTransformation::getA() {
    Rcpp::NumericMatrix A_curr(_AMatrix.get_nRow(), _AMatrix.get_nCol());

    for (int i = 0; i < A_curr.nrow(); ++i) {
        std::vector<double> curr_row = _AMatrix.get_Row(i);

        for (int j = 0; j < A_curr.ncol(); ++j) {
            A_curr(i, j) = curr_row[j];
        }
    }

    return A_curr;
}

Rcpp::NumericMatrix GibbsSamplerTransformation::getP() {
    Rcpp::NumericMatrix P_curr(_PMatrix.get_nRow(), _PMatrix.get_nCol());

    for (int i = 0; i < P_curr.nrow(); ++i) {
        std::vector<double> curr_row = _PMatrix.get_Row(i);

        for (int j = 0; j < P_curr.ncol(); ++j) {
            P_curr(i, j) = curr_row[j];
        }
    }

    return P_curr;
}

void GibbsSamplerTransformation::abc_mcmc(int burn, int iter, int thin, double tolerance) {
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

    //Rcpp::Rcout << "proposal\n";
    for (int s = 0; s < thin; ++s) {
        // propose theta
        _growth.propose(A_curr, P_curr);

        // check whether growth was accepted
        if (_growth.accepted) {
            Rcpp::Rcout << "Accepted growth proposal\n";
        }

        // update the A matrix
        Rcpp::Rcout << "weight: " << calcWeight() << "\n";
        weightAAtomicColumn(calcWeight());
        weightAColumn(calcWeight());

        // update the P matrix
        updatePRow();

    }
    //Rcpp::Rcout << "Done proposing\n";

    // save the theta's
    _theta(burn + iter - 1, Rcpp::_) = _growth.theta();
    _thresh(burn + iter - 1, 0) = _growth.rho;
    _thresh(burn + iter - 1, 1) = _growth.rho_thresh;
    
    // save the epsilon's
    _epsilon[burn + iter - 1] =_growth.epsilon();

    // track acceptance rate
    if (_growth.accepted == true) {
        accepted++;
    }

    proposed++;

}

// get relevant column and rows
Rcpp::NumericVector GibbsSamplerTransformation::getAcol() {
    unsigned int col = _AMatrix.get_nCol() - 1;
    std::vector<double> col_val = _AMatrix.get_Col(col);
    return Rcpp::wrap(col_val);
}

Rcpp::NumericVector GibbsSamplerTransformation::getProw() {
    unsigned int row = _PMatrix.get_nRow() - 1;
    std::vector<double> row_val = _PMatrix.get_Row(row);
    return Rcpp::wrap(row_val);
}

Rcpp::NumericVector GibbsSamplerTransformation::currentTheta() {
    return _growth.theta();
}
