#ifndef __COGAPS_INTERNAL_STATE_H__
#define __COGAPS_INTERNAL_STATE_H__

#include "Archive.h"
#include "Matrix.h"

#include <Rcpp.h>

typedef std::vector<Rcpp::NumericMatrix> SnapshotList;

enum GapsPhase
{
    GAPS_BURN,
    GAPS_COOL,
    GAPS_SAMP
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
    
    GapsInternalState(const Rcpp::NumericMatrix &D,
        const Rcpp::NumericMatrix &S, const Rcpp::NumericMatrix &fixedPatterns,
        const Rcpp::S4 &params)
            :
        chi2VecEquil(params.slot("nEquil")),
        nAtomsAEquil(params.slot("nEquil")),
        nAtomsPEquil(params.slot("nEquil")),
        chi2VecSample(params.slot("nSample")),
        nAtomsASample(params.slot("nSample")),
        nAtomsPSample(params.slot("nSample")),
        nIterA(10),
        nIterP(10),
        nEquil(params.slot("nEquil")),
        nEquilCool(params.slot("nEquilCool")),
        nSample(params.slot("nSample")),
        nSnapshots(params.slot("nSnapshots")),
        nOutputs(params.slot("nOutput")),
        messages(params.slot("messages")),
        iter(0),
        phase(GAPS_BURN),
        seed(params.slot("seed")), 
        checkpointInterval(params.slot("checkpointInterval")),
        numCheckpoints(0),
        nUpdatesA(0),
        nUpdatesP(0),
        sampler(D, S, params.slot("nFactor"), params.slot("alphaA"),
            params.slot("alphaP"), params.slot("maxGibbMassA"),
            params.slot("maxGibbMassP"), params.slot("singleCellRNASeq"),
            fixedPatterns, params.slot("whichMatrixFixed"))
    {}

    GapsInternalState(const Rcpp::NumericMatrix &D,
        const Rcpp::NumericMatrix &S, unsigned nF, unsigned nE, unsigned nS)
            :
        chi2VecEquil(nE), nAtomsAEquil(nE), nAtomsPEquil(nE),
        chi2VecSample(nS), nAtomsASample(nS), nAtomsPSample(nS),
        sampler(D, S, nF)
    {}
};

inline Archive& operator<<(Archive &ar, GapsInternalState &state)
{
    ar << state.chi2VecEquil << state.nAtomsAEquil << state.nAtomsPEquil
        << state.chi2VecSample << state.nAtomsASample << state.nAtomsPSample
        << state.nIterA << state.nIterP << state.nEquil << state.nEquilCool
        << state.nSample << state.nSnapshots << state.nOutputs << state.messages
        << state.iter << state.phase << state.seed << state.checkpointInterval
        << state.numCheckpoints << state.nUpdatesA << state.nUpdatesP
        << state.sampler;
    return ar;
}

inline Archive& operator>>(Archive &ar, GapsInternalState &state)
{
    ar >> state.chi2VecEquil >> state.nAtomsAEquil >> state.nAtomsPEquil
        >> state.chi2VecSample >> state.nAtomsASample >> state.nAtomsPSample
        >> state.nIterA >> state.nIterP >> state.nEquil >> state.nEquilCool
        >> state.nSample >> state.nSnapshots >> state.nOutputs >> state.messages
        >> state.iter >> state.phase >> state.seed >> state.checkpointInterval
        >> state.numCheckpoints >> state.nUpdatesA >> state.nUpdatesP
        >> state.sampler;
    return ar;
}

#endif