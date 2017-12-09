#ifndef __COGAPS_GIBBSSAMPLER_H__
#define __COGAPS_GIBBSSAMPLER_H__

#include "AtomicSupport.h"
#include "Matrix.h"

#include <vector>

class GibbsSampler
{
private:

    AtomicSupport mADomain, mPDomain;
    Matrix mDMatrix, mSMatrix, mAMatrix, mPMatrix, mAPMatrix;

    double mAnnealingTemp;
    double mChi2;

    bool mSingleCellRNASeq;

    bool death(char matrixLabel, const AtomicSupport &domain, const AtomicProposal &proposal);
    bool birth(char matrixLabel, const AtomicSupport &domain, const AtomicProposal &proposal);
    bool move(char matrixLabel, const AtomicSupport &domain, const AtomicProposal &proposal);
    bool exchange(char matrixLabel, const AtomicSupport &domain, const AtomicProposal &proposal);

    double logLikelihood();

    double computeDeltaLL(const MatrixChange &change);

    double getGibbsMass(char matrixLabel, double origMass, unsigned row, unsigned int col);

public:

    GibbsSampler(Rcpp::NumericMatrix D, Rcpp::NumericMatrix S, unsigned nFactor,
        double alphaA, double alphaP);

    void update(char matrixLabel);
};

#endif
