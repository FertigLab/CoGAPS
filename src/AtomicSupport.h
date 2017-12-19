#ifndef __COGAPS_ATOMIC_SUPPORT_H__
#define __COGAPS_ATOMIC_SUPPORT_H__

#include "Random.h"
#include "Matrix.h"

#include <map>
#include <fstream>
#include <vector>
#include <stdint.h>

struct AtomicProposal
{
    char label;
    char type;
    unsigned nChanges;

    uint64_t pos1;
    double delta1;

    uint64_t pos2;
    double delta2;

    AtomicProposal(char l, char t, uint64_t p1, double m1)
        : label(l), type(t), nChanges(1), pos1(p1), delta1(m1), pos2(0), delta2(0.0)
    {}

    AtomicProposal(char l, char t, uint64_t p1, double m1, uint64_t p2, double m2)
        : label(l), type(t), nChanges(2), pos1(p1), delta1(m1), pos2(p2), delta2(m2)
    {}
};

class AtomicSupport 
{
private:

    // label of this atomic domain (A/P)
    char mLabel;

    // storage of the atomic domain
    std::map<uint64_t, double> mAtomicDomain;
    uint64_t mNumAtoms;
    uint64_t mMaxNumAtoms;

    // total mass of all atoms
    double mTotalMass;

    // matrix information
    uint64_t mNumRows;
    uint64_t mNumCols;
    uint64_t mNumBins;
    uint64_t mBinSize;

    // average number of atoms per bin (must be > 0)
    double mAlpha;

    // expected magnitude of each atom (must be > 0)
    double mLambda;

#ifdef GAPS_DEBUG
    mutable std::vector<char> mProposalHistory;
    mutable std::vector<uint64_t> mAtomHistory;
#endif

#ifdef GAPS_INTERNAL_TESTS
public:
#endif

    // convert atomic position to row/col of the matrix
    uint64_t getRow(uint64_t pos) const;
    uint64_t getCol(uint64_t pos) const;

    // get atomic positions at random
    uint64_t randomFreePosition() const;
    uint64_t randomAtomPosition() const;

    // proposal helper functions
    AtomicProposal proposeBirth() const;
    AtomicProposal proposeDeath() const;
    AtomicProposal proposeMove() const;
    AtomicProposal proposeExchange() const;

    // update the mass of an atom, return the total amount changed
    double updateAtomMass(char type, uint64_t pos, double delta);

#ifndef GAPS_INTERNAL_TESTS
public:
#endif

    // constructor
    AtomicSupport(char label, uint64_t nrow, uint64_t ncol, double alpha=1.0,
        double lambda=1.0);

    // create and accept a proposal
    AtomicProposal makeProposal() const;
    MatrixChange acceptProposal(const AtomicProposal &prop);

    // convert an AtomicProposal to an ElementChange
    MatrixChange getMatrixChange(const AtomicProposal &prop) const;

    // getters
    double alpha() const {return mAlpha;}
    double lambda() const {return mLambda;}
    double totalMass() const {return mTotalMass;}
    uint64_t numAtoms() const {return mNumAtoms;}
    double at(uint64_t loc) const {return mAtomicDomain.at(loc);}

    // setters
    void setAlpha(double alpha) {mAlpha = alpha;}
    void setLambda(double lambda) {mLambda = lambda;}
    //void setMaxNumAtoms(uint64_t max) {mMaxNumAtoms = max;}

#ifdef GAPS_DEBUG
    std::vector<char> proposalHistory() {return mProposalHistory;}
    std::vector<uint64_t> atomHistory() {return mAtomHistory;}
#endif
};

#endif
