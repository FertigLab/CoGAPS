#include "GibbsSampler.h"

// clear all quantities related to the local matrix proposal
void GibbsSampler::clear_Proposal() {
    _Row_changed.clear();
    _Col_changed.clear();
    _mass_changed.clear();
    _atomicProposal.clear();
    _matrixElemChange.clear();
}

// clear all quantities related to the new local matrix proposal
void GibbsSampler::clear_new_Proposal() {
    _new_Row_changed.clear();
    _new_Col_changed.clear();
    _new_mass_changed.clear();
    _new_atomicProposal.clear();
    _new_matrixElemChange.clear();
}

void GibbsSampler::local_display_matrix2F(ofstream &outputFile, double **Mat_ptr,
        unsigned int n_row, unsigned int n_col) {
    for (unsigned int m = 0; m < n_row; ++m) {
        for (unsigned int n = 0; n < n_col; ++n) {
            outputFile << std::setw(10) << std::right;
            outputFile << Mat_ptr[m][n] << " ";
        }

        outputFile << endl;
    }
}

// -----------------------------------------------------------------------------
void GibbsSampler::output_atomicdomain(char atomic_label, unsigned long Samp_cycle) {
    char outputFilename[80];

    switch (atomic_label) {
        case 'A': {
            strcpy(outputFilename, _simulation_id.c_str());
            strcat(outputFilename, "_A_atomicdomain.txt");
            _AAtomicdomain.writeAtomicInfo(outputFilename, Samp_cycle);
            break;
        }

        case 'P': {
            strcpy(outputFilename, _simulation_id.c_str());
            strcat(outputFilename, "_P_atomicdomain.txt");
            _PAtomicdomain.writeAtomicInfo(outputFilename, Samp_cycle);
            break;
        }
    }
}

// -----------------------------------------------------------------------------
void GibbsSampler::compute_statistics_prepare_matrices(unsigned long statindx) {
    double **A = _AMatrix.get_matrix();
    double **P = _PMatrix.get_matrix();
    vector<double> k(_nFactor);  // normalized vector

    // compute the normalization vector
    for (int m = 0; m < _nFactor; ++m) {
        k[m] = 0.;

        for (int n = 0; n < _nCol; ++n) {
            k[m] += P[m][n];
        }

        if (k[m] == 0) { // when the whole row of P is zero, then don't do anything
            k[m] = 1.0;
        }
    }

    // construct the mean and var matrices at statindx = 1
    if (statindx == 1) {
        _Amean = new double * [_nRow];

        for (int m = 0; m < _nRow ; ++m) {
            _Amean[m] = new double [_nFactor];
        }

        for (int m = 0; m < _nRow; ++m) {
            for (int n = 0; n < _nFactor; ++n) {
                _Amean[m][n] = A[m][n] * k[n];
            }
        }

        _Asd = new double * [_nRow];

        for (int m = 0; m < _nRow ; ++m) {
            _Asd[m] = new double [_nFactor];
        }

        for (int m = 0; m < _nRow; ++m) {
            for (int n = 0; n < _nFactor; ++n) {
                _Asd[m][n] = pow(A[m][n] * k[n], 2);
            }
        }

        _Pmean = new double * [_nFactor];

        for (int m = 0; m < _nFactor ; ++m) {
            _Pmean[m] = new double [_nCol];
        }

        for (int m = 0; m < _nFactor; ++m) {
            for (int n = 0; n < _nCol; ++n) {
                _Pmean[m][n] = P[m][n] / k[m];
            }
        }

        _Psd = new double * [_nFactor];

        for (int m = 0; m < _nFactor ; ++m) {
            _Psd[m] = new double [_nCol];
        }

        for (int m = 0; m < _nFactor; ++m) {
            for (int n = 0; n < _nCol; ++n) {
                _Psd[m][n] = pow(P[m][n] / k[m] , 2);
            }
        }
    } // end of if-block for matrix construction statindx == 1

    // increment the mean and var matrices at statindx != 1
    if (statindx > 1) {
        for (int m = 0; m < _nRow; ++m) {
            for (int n = 0; n < _nFactor; ++n) {
                _Amean[m][n] += A[m][n] * k[n];
            }
        }

        for (int m = 0; m < _nRow; ++m) {
            for (int n = 0; n < _nFactor; ++n) {
                _Asd[m][n] += pow(A[m][n] * k[n], 2);
            }
        }

        for (int m = 0; m < _nFactor; ++m) {
            for (int n = 0; n < _nCol; ++n) {
                _Pmean[m][n] += P[m][n] / k[m];
            }
        }

        for (int m = 0; m < _nFactor; ++m) {
            for (int n = 0; n < _nCol; ++n) {
                _Psd[m][n] += pow(P[m][n] / k[m], 2);
            }
        }
    } // end of if-block for matrix incrementation statindx > 1
} // end of method compute_statistics_prepare_matrices

// -----------------------------------------------------------------------------
void GibbsSampler::compute_statistics(unsigned int Nstat,
                                      vector<vector <double> > &AMeanVect,
                                      vector<vector <double> > &AStdVect,
                                      vector<vector <double> > &PMeanVect,
                                      vector<vector <double> > &PStdVect
                                     ) {
    //Conor's Code for resizing the vectors corresponding to each matrix
    AMeanVect.resize(_nRow);

    for (int i = 0; i < _nRow; i++) {
        AMeanVect[i].resize(_nFactor);
    }

    AStdVect.resize(_nRow);

    for (int i = 0; i < _nRow; i++) {
        AStdVect[i].resize(_nFactor);
    }

    PMeanVect.resize(_nFactor);

    for (int i = 0; i < _nFactor; i++) {
        PMeanVect[i].resize(_nCol);
    }

    PStdVect.resize(_nFactor);

    for (int i = 0; i < _nFactor; i++) {
        PStdVect[i].resize(_nCol);
    }  //End vector resizing

    // compute statistics for A
    //Variable for holding temp Calculations
    double tempStat;

    //Changed to compute as vectors to pass into R
    for (int m = 0; m < _nRow ; ++m) {
        for (int n = 0; n < _nFactor; ++n) {
            tempStat = _Amean[m][n] / Nstat;
            _Amean[m][n] = tempStat;
            AMeanVect[m][n] = tempStat;
        }
    }

    for (int m = 0; m < _nRow ; ++m) {
        for (int n = 0; n < _nFactor; ++n) {
            tempStat = sqrt((_Asd[m][n] - Nstat * pow(_Amean[m][n], 2)) / (Nstat - 1));
            _Asd[m][n] = tempStat;
            AStdVect[m][n] = tempStat;
        }
    }

    // compute statistics for P
    for (int m = 0; m < _nFactor ; ++m) {
        for (int n = 0; n < _nCol; ++n) {
            tempStat = _Pmean[m][n] / Nstat;
            _Pmean[m][n] = tempStat;
            PMeanVect[m][n] = tempStat;
        }
    }

    for (int m = 0; m < _nFactor ; ++m) {
        for (int n = 0; n < _nCol; ++n) {
            tempStat = sqrt((_Psd[m][n] - Nstat * pow(_Pmean[m][n], 2)) / (Nstat - 1));
            _Psd[m][n] = tempStat;
            PStdVect[m][n] = tempStat;
        }
    }
}

// *****************************************************************************
// Adaptation from the original code:

bool GibbsSampler::checkOtherMatrix(char the_matrix_label, unsigned int iRow, unsigned int iCol,
                                    double const *const *otherMatrix) {
    unsigned int otherDim;

    // check that there is mass for this pattern in the P Matrix
    if (the_matrix_label == 'A') {
        for (otherDim = 0; otherDim < _nCol; otherDim++) {
            if (otherMatrix[iCol][otherDim] > epsilon) {
                return true;
            }
        }

        return false;
    }

    // check that there is mass for this pattern in the A Matrix
    else {
        for (otherDim = 0; otherDim < _nRow; otherDim++) {
            if (otherMatrix[otherDim][iRow] > epsilon) {
                return true;
            }
        }

        return false;
    }

    return false;
}

//-----------------------------------------------------------------

//-----------------------------------------------------------------
double GibbsSampler::getMass(char the_matrix_label, double origMass,
                             unsigned int iRow,
                             unsigned int iCol,
                             double const *const *otherMatrix,
                             double const *const *currentChainMatrix,
                             double const *const *D, double const *const *S,
                             double rng) {
    double DOUBLE_POSINF = numeric_limits<double>::max();
    unsigned int iGene = 0;
    unsigned int iPattern = 0;
    unsigned int iSample = 0;
    unsigned int jPattern = 0;
    double lambda;

    switch (the_matrix_label) {
        case 'A': {
            iGene = iRow;
            iPattern = iCol;
            lambda = _AAtomicdomain.getLambda();
            break;
        }

        case 'P': {
            iPattern = iRow;
            iSample = iCol;
            lambda = _PAtomicdomain.getLambda();
            break;
        }
    }

    // determine the parameters for finding the mass
    double s  = 0.;
    double su = 0.;
    double mock;

    switch (the_matrix_label) {
        case 'A': {
            //double Aeff;
            for (iSample = 0; iSample < _nCol; iSample++) {
                mock = D[iGene][iSample];

                for (jPattern = 0; jPattern < _nFactor; jPattern++) {
                    mock -= currentChainMatrix[iGene][jPattern] * otherMatrix[jPattern][iSample];
                }

                s += _annealingTemperature * pow(otherMatrix[iPattern][iSample], 2) /
                     (2 * pow(S[iGene][iSample], 2));
                su += _annealingTemperature * otherMatrix[iPattern][iSample] * mock /
                      (2 * pow(S[iGene][iSample], 2));
            }

            break;
        }

        case 'P': {
            for (iGene = 0; iGene < _nRow; iGene++) {
                mock = D[iGene][iSample];

                for (jPattern = 0; jPattern < _nFactor; jPattern++) {
                    mock -= otherMatrix[iGene][jPattern] * currentChainMatrix[jPattern][iSample];
                }

                s += _annealingTemperature * pow(otherMatrix[iGene][iPattern], 2) /
                     (2 * pow(S[iGene][iSample], 2));
                su += _annealingTemperature * otherMatrix[iGene][iPattern] * mock /
                      (2 * pow(S[iGene][iSample], 2));
            }

            break;
        }
    }

    double mean  = (2 * su - lambda) / (2 * s);
    double sd = 1. / sqrt(2 * s);
    // note: is bounded below by zero so have to use inverse sampling!
    // based upon algorithm in DistScalarRmath.cc (scalarRandomSample)
    double plower = sub_func::pnorm(0., mean, sd, DOUBLE_NEGINF, 0);
    double pupper = 1.;
    double u = plower + randgen('U', 0, 0) * (pupper - plower);
    // -------------------------------------------------------------
    // this line seems to be misplaced.
    // double newMass = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);
    double newMass = 0;
    // ------------------

    // if the likelihood is flat and nonzero,
    // force to sample strictly from the prior
    if (plower == 1 || s < 1.e-5 || newMass == DOUBLE_POSINF || pupper == 0) {
        if (origMass < 0) {    // death case
            newMass = abs(origMass);

        } else {
            newMass = 0.;  // birth case
        }
    } // end of first comparison

    else if (plower >= 0.99) {
        double tmp1 = sub_func::dnorm(0, mean, sd, false);
        double tmp2 = sub_func::dnorm(10 * lambda, mean, sd, false);

        if ((tmp1 > epsilon) && (fabs(tmp1 - tmp2) < epsilon))   {
            if (origMass < 0) {   // death case
                return 0.;
            }

            return origMass;   // birth case
        }
    } // end of second comparison

    else {
        newMass = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);  // both death and birth
    }  // end of if-block for the remaining case

    // limit the mass range
    switch (the_matrix_label) {
        case 'A': {
            if (newMass > _max_gibbsmassA) {
                newMass = _max_gibbsmassA;
            }

            break;
        }

        case 'P': {
            if (newMass > _max_gibbsmassP) {
                newMass = _max_gibbsmassP;
            }

            break;
        }
    }

    if (newMass < 0) {
        newMass = 0;    // due to the requirement that newMass > 0
    }

    return newMass;
}
