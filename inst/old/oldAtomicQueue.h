#ifndef __GAPS_ATOMIC_QUEUE_H
#define __GAPS_ATOMIC_QUEUE_H

#include <map>
#include <vector>

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

class AtomicQueue
{
private:

    // parameters
    char mLabel;
    float mAlpha;
    float mLambda;

    // storage of the atomic domain
    std::vector<Atom> mAtoms;
    std::map<uint64_t, uint64_t> mAtomPositions;
    uint64_t mNumAtoms;
    uint64_t mAtomCapacity;
    
    // matrix information
    unsigned mNumRows;
    unsigned mNumCols;
    unsigned mNumBins;
    unsigned mBinSize;

    // storage of proposal queue
    std::queue<AtomicProposal> mProposalQueue;
    std::vector<unsigned> mUsedIndices; // rows of A, cols of P
    unsigned mMinAtoms;
    unsigned mMaxAtoms;
    
    // proposal functions
    AtomicProposal proposeBirth() const;
    AtomicProposal proposeDeath() const;
    AtomicProposal proposeMove() const;
    AtomicProposal proposeExchange() const;

    // update the mass of an atom, return the total amount changed
    std::pair<float, bool> updateAtomMass(uint64_t pos, float delta);

    // functions for dealing with atomic data structure
    void addAtom(Atom atom); // O(logN)
    void removeAtom(uint64_t pos); // O(logN)
    AtomNeighbors getNeighbors(uint64_t pos) const; // O(logN)

    // get atomic positions at random
    uint64_t randomFreePosition() const; // O(1)
    uint64_t randomAtomPosition() const; // O(1)

public:

    // constructor
    AtomicQueue(char label, unsigned nrow, unsigned ncol);

    // parameter setters
    void setAlpha(float alpha) {mAlpha = alpha;}
    void setLambda(float lambda) {mLambda = lambda;}

    // queue operations
    void populate(unsigned limit);
    AtomicProposal pop_front();
    MatrixChange accept(const AtomicProposal &prop, MatrixChange &change);
    void reject(const AtomicProposal &prop);
    unsigned size() const;
    bool empty() const;

    // serialization
    friend Archive& operator<<(Archive &ar, AtomicQueue &sampler);
    friend Archive& operator>>(Archive &ar, AtomicQueue &sampler);
};

#endif