#ifndef __COGAPS_INTERNAL_STATE_H__
#define __COGAPS_INTERNAL_STATE_H__

#include "Archive.h"
#include "Matrix.h"

#include <Rcpp.h>

typedef std::vector<Rcpp::NumericMatrix> SnapshotList;

enum GapsPhase
{
    GAPS_CALIBRATION,
    GAPS_COOLING,
    GAPS_SAMPLING
};

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

    unsigned nSnapshots;
    unsigned nOutputs;
    bool messages;

    unsigned iter;
    GapsPhase phase;
    uint32_t seed;

    long checkpointInterval;
    unsigned numCheckpoints;

    unsigned nUpdatesA;
    unsigned nUpdatesP;

    GibbsSampler sampler;

    SnapshotList snapshotsA;
    SnapshotList snapshotsP;

    GapsInternalState(Rcpp::NumericMatrix DMatrix, Rcpp::NumericMatrix SMatrix,
        unsigned nFactor, float alphaA, float alphaP, unsigned nE,
        unsigned nEC, unsigned nS, float maxGibbsMassA,
        float maxGibbsMassP, Rcpp::NumericMatrix fixedPatterns,
        char whichMatrixFixed, bool msg, bool singleCellRNASeq,
        unsigned numOutputs, unsigned numSnapshots, uint32_t in_seed,
        unsigned cptInterval)
            :
        chi2VecEquil(nE), nAtomsAEquil(nE), nAtomsPEquil(nE),
        chi2VecSample(nS), nAtomsASample(nS), nAtomsPSample(nS),
        nIterA(10), nIterP(10), nEquil(nE), nEquilCool(nEC), nSample(nS),
        nSnapshots(numSnapshots), nOutputs(numOutputs), messages(msg),
        iter(0), phase(GAPS_CALIBRATION), seed(in_seed),
        checkpointInterval(cptInterval), numCheckpoints(0),
        nUpdatesA(0), nUpdatesP(0),
        sampler(DMatrix, SMatrix, nFactor, alphaA, alphaP,
            maxGibbsMassA, maxGibbsMassP, singleCellRNASeq, fixedPatterns,
            whichMatrixFixed)
    {}

    // empty internal state, ready to be loaded
    GapsInternalState(unsigned nE, unsigned nS, unsigned nRow, unsigned nCol,
    unsigned nFactor)
            :
        chi2VecEquil(nE), nAtomsAEquil(nE), nAtomsPEquil(nE),
        chi2VecSample(nS), nAtomsASample(nS), nAtomsPSample(nS),
        sampler(nRow, nCol, nFactor)
    {}
};

inline void operator<<(Archive &ar, GapsInternalState &state)
{
    ar << state.chi2VecEquil;
    ar << state.nAtomsAEquil;
    ar << state.nAtomsPEquil;
    ar << state.chi2VecSample;
    ar << state.nAtomsASample;
    ar << state.nAtomsPSample;
    ar << state.nIterA;
    ar << state.nIterP;
    ar << state.nEquil;
    ar << state.nEquilCool;
    ar << state.nSample;
    ar << state.nSnapshots;
    ar << state.nOutputs;
    ar << state.messages;
    ar << state.iter;
    ar << state.phase;
    ar << state.seed;
    ar << state.checkpointInterval;
    ar << state.numCheckpoints;
    ar << state.nUpdatesA;
    ar << state.nUpdatesP;
    ar << state.sampler;
}

inline void operator>>(Archive &ar, GapsInternalState &state)
{
    ar >> state.chi2VecEquil;
    ar >> state.nAtomsAEquil;
    ar >> state.nAtomsPEquil;
    ar >> state.chi2VecSample;
    ar >> state.nAtomsASample;
    ar >> state.nAtomsPSample;
    ar >> state.nIterA;
    ar >> state.nIterP;
    ar >> state.nEquil;
    ar >> state.nEquilCool;
    ar >> state.nSample;
    ar >> state.nSnapshots;
    ar >> state.nOutputs;
    ar >> state.messages;
    ar >> state.iter;
    ar >> state.phase;
    ar >> state.seed;
    ar >> state.checkpointInterval;
    ar >> state.numCheckpoints;
    ar >> state.nUpdatesA;
    ar >> state.nUpdatesP;
    ar >> state.sampler;
}

#endif