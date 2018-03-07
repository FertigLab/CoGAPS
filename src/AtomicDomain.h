#ifndef __GAPS_ATOMIC_DOMAIN_H__
#define __GAPS_ATOMIC_DOMAIN_H__

#include "Archive.h"

#include <boost/unordered_set.hpp>
#include <stdint.h>
#include <cstddef>
#include <vector>
#include <map>

struct Atom
{
    uint64_t pos;
    float mass;

    Atom* left;
    Atom* right;
    
    Atom(uint64_t p, float m)
        : pos(p), mass(m), left(NULL), right(NULL)
    {}

    bool operator==(const Atom &other) const
    {
        return pos == other.pos;
    }
};

// data structure that holds atoms
class AtomicDomain
{
private:

    // domain storage
    std::vector<Atom> mAtoms;
    std::map<uint64_t, uint64_t> mAtomPositions;

    // TODO google_dense_set - first profile and benchmark
    boost::unordered_set<uint64_t> mUsedPositions;

public:

    // access atoms
    Atom front() const;
    Atom randomAtom() const;
    uint64_t randomFreePosition() const;
    unsigned size() const;

    // modify domain
    void insert(uint64_t pos, float mass);
    void erase(uint64_t pos);
    void updateMass(uint64_t pos, float newMass);

    // serialization
    friend Archive& operator<<(Archive &ar, AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);
};

#endif
