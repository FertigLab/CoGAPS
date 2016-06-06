#include "GibbsSampler.h"

// -----------------------------------------------------------------------------
unsigned long long atomicSize = std::numeric_limits<unsigned long long>::max();
// -----------------------------------------------------------------------------

// ******************** CONSTRUCTOR ********************************************
GibbsSampler:: GibbsSampler(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                            double alphaA, double alphaP, double nMaxA, double nMaxP,
                            unsigned long nIterA, unsigned long nIterP,
                            double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                            unsigned long long atomicSize,
                            char label_A, char label_P, char label_D, char label_S,
                            vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                            const string &simulation_id)
    : _DMatrix(DVector, label_D),
      _SMatrix(SVector, label_S) {
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


// *************** METHODS FOR INITIALIZATION, DISPLAY, OUTPUT ***********************
void GibbsSampler::init_AMatrix_and_PMatrix() {
    // extract information from D as parameters
    _nRow = _DMatrix.get_nRow();
    _nCol = _DMatrix.get_nCol();
    // initialize matrices A and p
    _AMatrix.born_matrix(_nRow, _nFactor, _label_A, _alphaA);
    _PMatrix.born_matrix(_nFactor, _nCol, _label_P, _alphaP);
}

void GibbsSampler::init_AAtomicdomain_and_PAtomicdomain() {
    // extract information from D as parameters
    _nRow = _DMatrix.get_nRow();
    _nCol = _DMatrix.get_nCol();
    double D_mean = _DMatrix.cal_mean();
    // calcuate #Bins and lambda for the atomic spaces
    _nBinsA = _nRow * _nFactor;
    _lambdaA = _alphaA * sqrt(_nFactor / D_mean) * _lambdaA_scale_factor;
    _nBinsP = _nFactor * _nCol;
    _lambdaP = _alphaP * sqrt(_nFactor / D_mean) * _lambdaP_scale_factor;
    // calculate the maximum gibbs mass for A and p
    _max_gibbsmassA = _max_gibbsmass_paraA / _lambdaA;
    _max_gibbsmassP = _max_gibbsmass_paraP / _lambdaP;
    // initialize the atomic spaces
    _AAtomicdomain.initializeAtomic(_nBinsA, atomicSize, _alphaA, _lambdaA, _label_A);
    _PAtomicdomain.initializeAtomic(_nBinsP, atomicSize, _alphaP, _lambdaP, _label_P);
    //cout << "_lambdaA = " << _lambdaA << ", _max_gibbsmassA = " << _max_gibbsmassA << endl;
    //cout << "_lambdaP = " << _lambdaP << ", _max_gibbsmassP = " << _max_gibbsmassP << endl << endl;
}

// For fixing one domain in R
void GibbsSampler::init_AAtomicdomain_and_PAtomicdomain(char fixeddomain, vector<vector<double> > ReadBinProbs) {
    // extract information from D as parameters
    _nRow = _DMatrix.get_nRow();
    _nCol = _DMatrix.get_nCol();
    double D_mean = _DMatrix.cal_mean();
    // calcuate #Bins and lambda for the atomic spaces
    _nBinsA = _nRow * _nFactor;
    _lambdaA = _alphaA * sqrt(_nFactor / D_mean) * _lambdaA_scale_factor;
    _nBinsP = _nFactor * _nCol;
    _lambdaP = _alphaP * sqrt(_nFactor / D_mean) * _lambdaP_scale_factor;
    // calculate the maximum gibbs mass for A and p
    _max_gibbsmassA = _max_gibbsmass_paraA / _lambdaA;
    _max_gibbsmassP = _max_gibbsmass_paraP / _lambdaP;

    // initialize the atomic spaces (fixed or not)
    if (fixeddomain == 'A') {
        _AAtomicdomain.FixedBins_initializeAtomic(_nBinsA, atomicSize, _alphaA, _lambdaA, _label_A, ReadBinProbs);
        _PAtomicdomain.initializeAtomic(_nBinsP, atomicSize, _alphaP, _lambdaP, _label_P);

    } else {
        _AAtomicdomain.initializeAtomic(_nBinsA, atomicSize, _alphaA, _lambdaA, _label_A);
        _PAtomicdomain.FixedBins_initializeAtomic(_nBinsP, atomicSize, _alphaP, _lambdaP, _label_P, ReadBinProbs);
    }
}

// For fixing two domains in R
void GibbsSampler::init_AAtomicdomain_and_PAtomicdomain(vector<vector<double> > ReadBinProbsA,
        vector<vector<double> > ReadBinProbsP) {
    // extract information from D as parameters
    _nRow = _DMatrix.get_nRow();
    _nCol = _DMatrix.get_nCol();
    double D_mean = _DMatrix.cal_mean();
    // calcuate #Bins and lambda for the atomic spaces
    _nBinsA = _nRow * _nFactor;
    _lambdaA = _alphaA * sqrt(_nFactor / D_mean) * _lambdaA_scale_factor;
    _nBinsP = _nFactor * _nCol;
    _lambdaP = _alphaP * sqrt(_nFactor / D_mean) * _lambdaP_scale_factor;
    // calculate the maximum gibbs mass for A and p
    _max_gibbsmassA = _max_gibbsmass_paraA / _lambdaA;
    _max_gibbsmassP = _max_gibbsmass_paraP / _lambdaP;
    // initialize the atomic spaces (BOTH FIXED)
    _AAtomicdomain.FixedBins_initializeAtomic(_nBinsA, atomicSize, _alphaA, _lambdaA, _label_A, ReadBinProbsA);
    _PAtomicdomain.FixedBins_initializeAtomic(_nBinsP, atomicSize, _alphaP, _lambdaP, _label_P, ReadBinProbsP);
}
