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
    float delta1;

    uint64_t pos2;
    float delta2;

    AtomicProposal(char l, char t, uint64_t p1, float m1)
        : label(l), type(t), nChanges(1), pos1(p1), delta1(m1), pos2(0), delta2(0.f)
    {}

    AtomicProposal(char l, char t, uint64_t p1, float m1, uint64_t p2, float m2)
        : label(l), type(t), nChanges(2), pos1(p1), delta1(m1), pos2(p2), delta2(m2)
    {}
};

struct Atom
{
    uint64_t pos;
    float mass;
    
    Atom(uint64_t p, float m) : pos(p), mass(m) {}
};

struct AtomNeighbors
{
    Atom left;
    Atom right;

    AtomNeighbors(uint64_t lp, float lm, uint64_t rp, float rm)
        : left(Atom(lp, lm)), right(Atom(rp, rm))
    {}
};

class AtomicSupport 
{
private:

#ifdef GAPS_INTERNAL_TESTS
public:
#endif

    // label of this atomic domain (A/P)
    char mLabel;

    // storage of the atomic domain
    std::vector<Atom> mAtoms;
    std::map<uint64_t, uint64_t> mAtomPositions;
    uint64_t mNumAtoms;
    uint64_t mMaxNumAtoms;

    // total mass of all atoms
    float mTotalMass;

    // matrix information
    uint64_t mNumRows;
    uint64_t mNumCols;
    uint64_t mNumBins;
    uint64_t mBinSize;

    // average number of atoms per bin (must be > 0)
    float mAlpha;

    // expected magnitude of each atom (must be > 0)
    float mLambda;

    // functions for dealing with atomic data structure
    void addAtom(Atom atom); // O(logN)
    void removeAtom(uint64_t pos); // O(logN)
    AtomNeighbors getNeighbors(uint64_t pos) const; // O(logN)

    // get atomic positions at random
    uint64_t randomFreePosition() const;
    uint64_t randomAtomPosition() const;

    // proposal helper functions
    AtomicProposal proposeBirth() const;
    AtomicProposal proposeDeath() const;
    AtomicProposal proposeMove() const;
    AtomicProposal proposeExchange() const;

    // update the mass of an atom, return the total amount changed
    float updateAtomMass(uint64_t pos, float delta);

public:

    // constructor
    AtomicSupport(char label, uint64_t nrow, uint64_t ncol, float alpha=1.0,
        float lambda=1.0);

    // create and accept a proposal
    AtomicProposal makeProposal() const;
    MatrixChange acceptProposal(const AtomicProposal &prop);

    // convert an AtomicProposal to an ElementChange
    MatrixChange getMatrixChange(const AtomicProposal &prop) const;

    // convert atomic position to row/col of the matrix
    uint64_t getRow(uint64_t pos) const;
    uint64_t getCol(uint64_t pos) const;

    // getters
    float alpha() const {return mAlpha;}
    float lambda() const {return mLambda;}
    float totalMass() const {return mTotalMass;}
    uint64_t numAtoms() const {return mNumAtoms;}
    float at(uint64_t pos) const {return mAtoms[mAtomPositions.at(pos)].mass;}

    // setters
    void setAlpha(float alpha) {mAlpha = alpha;}
    void setLambda(float lambda) {mLambda = lambda;}

    friend Archive& operator<<(Archive &ar, AtomicSupport &sampler);
    friend Archive& operator>>(Archive &ar, AtomicSupport &sampler);
};

#endif
