#ifndef __COGAPS_INTERNAL_STATE_H__
#define __COGAPS_INTERNAL_STATE_H__

#include "Archive.h"
#include "Matrix.h"
#include "GibbsSampler.h"

#include <Rcpp.h>

typedef std::vector<Rcpp::NumericMatrix> SnapshotList;

enum GapsPhase
{
    GAPS_BURN,
    GAPS_COOL,
    GAPS_SAMP
};

// maybe cleaner if this was a class with member functions? GapsRunner
struct GapsInternalState
{
    Vector chi2VecEquil;
    Vector nAtomsAEquil;
    Vector nAtomsPEquil;

    Vector chi2VecSample;
    Vector nAtomsASample;
    Vector nAtomsPSample;

    unsigned nIterA;
    unsigned nIterP;
    
    unsigned nEquil;
    unsigned nEquilCool;
    unsigned nSample;
    unsigned nFactor;

    unsigned nSnapshots;
    unsigned nOutputs;
    bool messages;

    unsigned iter;
    GapsPhase phase;
    uint32_t seed;

    long checkpointInterval;

    unsigned nUpdatesA;
    unsigned nUpdatesP;
    
    unsigned nPumpSamples;

    AmplitudeGibbsSampler ASampler;
    PatternGibbsSampler PSampler;
    
    SnapshotList snapshotsA;
    SnapshotList snapshotsP;

    GapsInternalState(const Rcpp::NumericMatrix &D,
        const Rcpp::NumericMatrix &S, unsigned nF, unsigned nE, unsigned nEC,
        unsigned nS, unsigned nOut, unsigned nSnap, float alphaA, float alphaP,
        float maxGibbsMassA, float maxGibbsMassP, int sd, bool msgs,
        bool singleCellRNASeq, char whichMatrixFixed,
        const Rcpp::NumericMatrix &FP, unsigned cptInterval)
        //PumpThreshold pumpThreshold, unsigned numPumpSamples)
            :
        chi2VecEquil(nE), nAtomsAEquil(nE), nAtomsPEquil(nE),
        chi2VecSample(nS), nAtomsASample(nS), nAtomsPSample(nS),
        nIterA(10), nIterP(10), nEquil(nE), nEquilCool(nEC), nSample(nS),
        nSnapshots(nSnap), nOutputs(nOut), messages(msgs), iter(0),
        phase(GAPS_BURN), seed(sd), checkpointInterval(cptInterval),
        nUpdatesA(0), nUpdatesP(0), //nPumpSamples(numPumpSamples),
        ASampler(D, S, nF, alphaA, maxGibbsMassA),
        PSampler(D, S, nF, alphaP, maxGibbsMassP)
    {}

    GapsInternalState(const Rcpp::NumericMatrix &D,
        const Rcpp::NumericMatrix &S, unsigned nF, unsigned nE, unsigned nS)
            :
        chi2VecEquil(nE), nAtomsAEquil(nE), nAtomsPEquil(nE),
        chi2VecSample(nS), nAtomsASample(nS), nAtomsPSample(nS),
        ASampler(D, S, nF), PSampler(D, S, nF)
    {}
};

inline Archive& operator<<(Archive &ar, GapsInternalState &state)
{
/*
    ar << state.chi2VecEquil << state.nAtomsAEquil << state.nAtomsPEquil
        << state.chi2VecSample << state.nAtomsASample << state.nAtomsPSample
        << state.nIterA << state.nIterP << state.nEquil << state.nEquilCool
        << state.nSample << state.nSnapshots << state.nOutputs << state.messages
        << state.iter << state.phase << state.seed << state.checkpointInterval
        << state.nUpdatesA << state.nUpdatesP << state.nPumpSamples;
        //<< state.sampler;
*/
    return ar;
}

inline Archive& operator>>(Archive &ar, GapsInternalState &state)
{
/*  
    ar >> state.chi2VecEquil >> state.nAtomsAEquil >> state.nAtomsPEquil
        >> state.chi2VecSample >> state.nAtomsASample >> state.nAtomsPSample
        >> state.nIterA >> state.nIterP >> state.nEquil >> state.nEquilCool
        >> state.nSample >> state.nSnapshots >> state.nOutputs >> state.messages
        >> state.iter >> state.phase >> state.seed >> state.checkpointInterval
        >> state.nUpdatesA >> state.nUpdatesP >> state.nPumpSamples;
        //>> state.sampler;
*/
    return ar;
}

#endif