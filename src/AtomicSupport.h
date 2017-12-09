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
    char type;
    unsigned nChanges;

    uint64_t pos1;
    double delta1;

    uint64_t pos2;
    double delta2;

    AtomicProposal(char t, uint64_t p1, double m1)
        : type(t), nChanges(1), pos1(p1), delta1(m1), pos2(0), delta2(0.0)
    {}

    AtomicProposal(char t, uint64_t p1, double m1, uint64_t p2, double m2)
        : type(t), nChanges(2), pos1(p1), delta1(m1), pos2(p2), delta2(m2)
    {}
};

class AtomicSupport 
{
private:

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

    // average number of atoms per bin
    double mAlpha;

    // expected magnitude of each atom
    double mLambda;     

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

    // update the mass of an atom
    void updateAtomMass(uint64_t pos, double delta);

public:

    // constructor
    AtomicSupport(uint64_t nrow, uint64_t ncol);

    // create and accept a proposal
    AtomicProposal makeProposal() const;
    void acceptProposal(const AtomicProposal &prop);

    // convert an AtomicProposal to an ElementChange
    MatrixChange getMatrixChange(const AtomicProposal &prop) const;

    // write atomic domain to file
    void write(const std::string &outputFilename, bool append) const;
    
    // getters
    double alpha() const {return mAlpha;}
    double lambda() const {return mLambda;}
    double totalMass() const {return mTotalMass;}
    uint64_t numAtoms() const {return mNumAtoms;}

    // setters
    void setAlpha(double alpha) {mAlpha = alpha;}
    void setLambda(double lambda) {mLambda = lambda;}
    void setMaxNumAtoms(uint64_t max) {mMaxNumAtoms = max;}

    // TODO remove support for this
    std::map<uint64_t, double> getDomain() {return mAtomicDomain;}
};

#endif
