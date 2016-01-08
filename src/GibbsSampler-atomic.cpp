#include "GibbsSampler.h"

// ********* METHODS TO GO BETWEEN ATOMIC SPACE AND MATRIX  ********************
unsigned int GibbsSampler::getRow(char matrix_label , unsigned int iBin) {
    switch (matrix_label) {
        case 'A': { // A - horizontal addressing implicit
            return (floor(iBin / _nFactor));
        }

        case 'P': { // P - vertical addressing implicit
            return iBin % _nFactor;
        }
    }

    // EJF dummy return to avoid warnings
    return iBin;
}

unsigned int GibbsSampler::getCol(char matrix_label , unsigned int iBin) {
    switch (matrix_label) {
        case 'A': { // A - horizontal addressing implicit
            return iBin % _nFactor;
        }

        case 'P': { // P - vertical addressing implicit
            return floor(iBin / _nFactor);
        }
    }

    // EJF dummy return to avoid warnings
    return iBin;
}


unsigned int GibbsSampler::getTotNumAtoms(char matrix_label) {
    switch (matrix_label) {
        case 'A': {
            return _AAtomicdomain.getNAtom();
        }

        case 'P': {
            return _PAtomicdomain.getNAtom();
        }
    }

    // EJF dummy return to avoid warnings
    return 0;
}

vector <vector <vector <double> > > GibbsSampler::getNormedMatrices() {
    double **A = _AMatrix.get_matrix();
    double **P = _PMatrix.get_matrix();
    vector <vector <double> > AMatrixNormed;
    AMatrixNormed.resize(_nRow, vector<double>(_nFactor, 0.0));
    vector <vector <double> > PMatrixNormed;
    PMatrixNormed.resize(_nFactor, vector<double>(_nCol, 0.0));
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

    for (int m = 0; m < _nRow; ++m) {
        for (int n = 0; n < _nFactor; ++n) {
            AMatrixNormed[m][n] = A[m][n] * k[n];
        }
    }

    for (int m = 0; m < _nFactor; ++m) {
        for (int n = 0; n < _nCol; ++n) {
            PMatrixNormed[m][n] = P[m][n] / k[m];
        }
    }

    vector <vector <vector <double> > > NormedMatrices;
    NormedMatrices.push_back(AMatrixNormed);
    NormedMatrices.push_back(PMatrixNormed);
    return NormedMatrices;
}

// ======CoGAPS Test methods============

// Get the full atomic domain from the atomic space
map <unsigned long long, double> GibbsSampler::getAtomicDomain(char matrix_label) {
    map <unsigned long long, double> zero;

    if (matrix_label == 'A') {
        return _AAtomicdomain.getDomain();

    } else if (matrix_label == 'P') {
        return _PAtomicdomain.getDomain();
    }

    return zero;
}

// Manually calculate the matrix A from the atomic space passed in.
vector <vector <double> > GibbsSampler::createSampleAMat(map <unsigned long long, double> ADomain) {
    // Make that a matrix
    vector <vector <double> > SampleAMatrix;
    SampleAMatrix.resize(_nRow);

    for (int i = 0; i < _nRow; i++) {
        SampleAMatrix[i].resize(_nFactor, 0);
    }

    //Copy the parsed domain
    map<unsigned long long, double>::const_iterator iter;

    for (iter = ADomain.begin(); iter != ADomain.end(); iter++) {
        unsigned int theBin = _AAtomicdomain.getBin(iter->first);
        // Put the mass in the bin in the matrix
        unsigned int theRow = getRow('A', theBin);
        unsigned int theCol = getCol('A', theBin);
        SampleAMatrix[theRow][theCol] += iter->second;
    }

    return SampleAMatrix;
}

// Manually calculate the matrix P from the atomic space passed in.
vector <vector <double> > GibbsSampler::createSamplePMat(map <unsigned long long, double> PDomain) {
    // Make that a matrix
    vector <vector <double> > SamplePMatrix;
    SamplePMatrix.resize(_nFactor);

    for (int i = 0; i < _nFactor; i++) {
        SamplePMatrix[i].resize(_nCol, 0);
    }

    //Copy the parsed domain
    map<unsigned long long, double>::const_iterator iter;

    for (iter = PDomain.begin(); iter != PDomain.end(); iter++) {
        unsigned int theBin = _PAtomicdomain.getBin(iter->first);
        // Put the mass in the bin in the matrix
        unsigned int theRow = getRow('P', theBin);
        unsigned int theCol = getCol('P', theBin);
        SamplePMatrix[theRow][theCol] += iter->second;
    }

    return SamplePMatrix;
}

// Manually calculate the chi squared value based on the 2 matrices passed in
double GibbsSampler::ManualCalcChiSqu(vector <vector <double> > SampleAMat, vector <vector <double> > SamplePMat) {
    Matrix SampleAMatrix(SampleAMat, 'S');
    Matrix SamplePMatrix(SamplePMat, 'S');
    double **D = _DMatrix.get_matrix();
    double **S = _SMatrix.get_matrix();
    double **A = SampleAMatrix.get_matrix();
    double **P = SamplePMatrix.get_matrix();
    return GAPSNorm::calChi2(D, S, A, P, _nRow, _nCol, _nFactor);
}

