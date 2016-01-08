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
void GibbsSampler::check_resultsF(ofstream &outputFile) {
    //EJF double const * const * D = _DMatrix.get_matrix();
    //EJF double const * const * S = _SMatrix.get_matrix();
    double const *const *A = _AMatrix.get_matrix();
    double const *const *P = _PMatrix.get_matrix();
    vector<vector<double> > AP;
    AP.resize(_nRow, vector<double>(_nCol, 0.0));

    for (unsigned int m = 0; m < _nRow; ++m) {
        for (unsigned int n = 0; n < _nCol; ++n) {
            for (unsigned int k = 0; k < _nFactor; ++k) {
                AP[m][n] += A[m][k] * P[k][n];
            }
        }
    }

    outputFile << "The product matrix AP = A*P is: " << endl;
    outputFile << endl;

    for (unsigned int m = 0; m < _nRow; ++m) {
        for (unsigned int n = 0; n < _nCol; ++n) {
            outputFile << setiosflags(ios::right) << setw(10) << AP[m][n] << " ";
        }

        outputFile << endl;
    }

    outputFile << endl;
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
