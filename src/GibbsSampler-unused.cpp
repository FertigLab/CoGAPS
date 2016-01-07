#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "GibbsSampler.h"

using namespace std;
using namespace gaps;
using std::vector;

// -----------------------------------------------------------------------------
unsigned long long atomicSize = std::numeric_limits<unsigned long long>::max();
const double DOUBLE_POSINF = std::numeric_limits<double>::max();
const double DOUBLE_NEGINF = -std::numeric_limits<double>::max();
const double epsilon = 1e-10;
// -----------------------------------------------------------------------------

GibbsSampler:: GibbsSampler(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                            double alphaA, double alphaP, double nMaxA, double nMaxP,
                            unsigned long nIterA, unsigned long nIterP,
                            double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                            unsigned long long atomicSize,
                            char label_A, char label_P, char label_D, char label_S,
                            const string &datafile, const string &variancefile,
                            const string &simulation_id)
    : _DMatrix(datafile.c_str(), label_D),
      _SMatrix(variancefile.c_str(), label_S) {
    _nEquil = nEquil;
    _nSample = nSample;
    _nFactor = nFactor;
    _alphaA = alphaA;
    _alphaP = alphaP;
    _nMaxA = nMaxA;
    _nMaxP = nMaxP;
    _nIterA = nIterA;
    _nIterP = nIterP;
    _max_gibbsmass_paraA = max_gibbsmass_paraA;
    _max_gibbsmass_paraP = max_gibbsmass_paraP;
    _lambdaA_scale_factor = 1;
    _lambdaP_scale_factor = 1;
    _simulation_id = simulation_id;
    _atomicSize = atomicSize;
    _label_A = label_A;
    _label_P = label_P;
    _iter = 1;  // tmp use
    _annealingTemperature = 0.4; // tmp use
    _sysChi2 = 0.0; // tmp use
}

// -----------------------------------------------------------------------------
void GibbsSampler::output_computing_info(char outputFilename[],
        unsigned long Equil_cycle, unsigned long nEquil,
        unsigned long Samp_cycle, unsigned long nSample) {
    ofstream outputFile;
    outputFile.open(outputFilename, ios::out | ios::app);
    outputFile << " *************************************************** " << endl;
    outputFile << " --------------- NEW ROUND ------------------------- " << endl;
    outputFile << " *************************************************** " << endl << endl;
    outputFile << "Equilibration cycle index = " << Equil_cycle << endl;
    outputFile << "Total number of equilibrating cycles to perform = " <<  nEquil << endl;
    outputFile << "Sampling cycle index = " << Samp_cycle << endl;
    outputFile << "Total number of sampling cycles to perform = " <<  nSample << endl;
    outputFile << "System Chi2-value = " << _sysChi2 << endl;
    _AMatrix.display_matrixF(outputFile);
    _PMatrix.display_matrixF(outputFile);
    check_resultsF(outputFile);
    outputFile.close();
}

void GibbsSampler::cal_delloglikelihood_example() {
}

// ---------------- For checking against computeDeltaLL
double GibbsSampler::computeDeltaLL2(char the_matrix_label,
                                     double const *const *D,
                                     double const *const *S,
                                     double const *const *A,
                                     double const *const *P,
                                     unsigned int the_nChange_matrixElemChange,
                                     const vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange) {
    double DelLL = 0.0;

    switch (the_matrix_label) {
        case 'A': {
            DelLL = GAPSNorm::calcDeltaLLGen('A', D, S, A, P, the_matrixElemChange, _nRow,
                                             _nCol, _nFactor);
            break;
        }

        case 'P': {
            DelLL = GAPSNorm::calcDeltaLLGen('P', D, S, A, P, the_matrixElemChange, _nRow,
                                             _nCol, _nFactor);
            break;
        }
    }

    return DelLL;
} // end of computeDeltaLL2
// -------------------------------------

void GibbsSampler::update_example(char atomic_domain_label) {
} // end of update_example

// *************** METHODS FOR MAKING PROPOSAL ********************************
// -----------------------------------------------------------------------------
vector<vector<double> > GibbsSampler::atomicProposal2Matrix(char atomic_domain_label,
        double const *const *origMatrix) {
    unsigned int bin;
    unsigned int chRow, chCol;

    switch (atomic_domain_label) {
        case 'A': {
            vector<vector<double> > newMatrix(_nRow, vector<double>(_nFactor, 0));
            map<unsigned long long, double> proposal = _AAtomicdomain.getProposedAtoms();

            for (map<unsigned long long, double>::const_iterator
                    iter = proposal.begin(); iter != proposal.end(); ++ iter) {
                bin = _AAtomicdomain.getBin(iter->first);
                chRow = getRow('A', bin);
                chCol = getCol('A', bin);
                newMatrix[chRow][chCol] += iter->second;
            }

            return newMatrix;
        } // end of case 'A'

        case 'P': {
            vector<vector<double> > newMatrix(_nFactor, vector<double>(_nCol, 0));
            map<unsigned long long, double> proposal = _PAtomicdomain.getProposedAtoms();

            for (map<unsigned long long, double>::const_iterator
                    iter = proposal.begin(); iter != proposal.end(); ++ iter) {
                bin = _PAtomicdomain.getBin(iter->first);
                chRow = getRow('P', bin);
                chCol = getCol('P', bin);
                newMatrix[chRow][chCol] += iter->second;
            }

            return newMatrix;
        } // end of case 'P'
    } // end of switch

    // EJF dummy return to avoid warning
    vector<vector<double> > newMatrix(0, vector<double>(0, 0));
    return newMatrix;
} // end of atomicProposal2Matrix

// -----------------------------------------------------------------------------
vector<vector<double> > GibbsSampler::atomicProposal2FullMatrix(char atomic_domain_label,
        double const *const *origMatrix) {
    unsigned int bin;
    unsigned int chRow, chCol;

    switch (atomic_domain_label) {
        case 'A': {
            vector<vector<double> > FullnewMatrix(_nRow, vector<double>(_nFactor, 0));

            for (unsigned int iRow = 0; iRow < _nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < _nFactor; ++iCol) {
                    if (origMatrix[iRow][iCol] < epsilon) {
                        FullnewMatrix[iRow][iCol] = 0.;

                    } else {
                        FullnewMatrix[iRow][iCol] = origMatrix[iRow][iCol];
                    }
                }
            } // end of for block

            map<unsigned long long, double> proposal = _AAtomicdomain.getProposedAtoms();

            for (map<unsigned long long, double>::const_iterator
                    iter = proposal.begin(); iter != proposal.end(); ++ iter) {
                bin = _AAtomicdomain.getBin(iter->first);
                chRow = getRow('A', bin);
                chCol = getCol('A', bin);
                FullnewMatrix[chRow][chCol] += iter->second;
            }

            return FullnewMatrix;
        } // end of case 'A'

        case 'P': {
            vector<vector<double> > FullnewMatrix(_nFactor, vector<double>(_nCol, 0));

            for (unsigned int iRow = 0; iRow < _nFactor; ++iRow) {
                for (unsigned int iCol = 0; iCol < _nCol; ++iCol) {
                    if (origMatrix[iRow][iCol] < epsilon) {
                        FullnewMatrix[iRow][iCol] = 0.;

                    } else {
                        FullnewMatrix[iRow][iCol] = origMatrix[iRow][iCol];
                    }
                }
            } // end of for block

            map<unsigned long long, double> proposal = _PAtomicdomain.getProposedAtoms();

            for (map<unsigned long long, double>::const_iterator
                    iter = proposal.begin(); iter != proposal.end(); ++ iter) {
                bin = _PAtomicdomain.getBin(iter->first);
                chRow = getRow('P', bin);
                chCol = getCol('P', bin);
                FullnewMatrix[chRow][chCol] += iter->second;
            }

            return FullnewMatrix;
        } // end of case 'P'
    } // end of switch

    // EJF dummy return to avoid warning
    vector<vector<double> > newMatrix(0, vector<double>(0, 0));
    return newMatrix;
} // end of atomicProposal2FullMatrix

double GibbsSampler::get_AnnealingTemperature() {
    return _annealingTemperature;
}

// -----------------------------------------------------------------------------
void GibbsSampler::detail_check(char outputchi2_Filename[]) {
    double chi2 = 2.*cal_logLikelihood();
    ofstream outputchi2_File;
    outputchi2_File.open(outputchi2_Filename, ios::out | ios::app);
    outputchi2_File << chi2 << endl;
    outputchi2_File.close();
}
