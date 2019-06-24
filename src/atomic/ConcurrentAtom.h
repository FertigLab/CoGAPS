#ifndef __COGAPS_CONCURRENT_ATOM_H__
#define __COGAPS_CONCURRENT_ATOM_H__

struct ConcurrentAtom;
class ConcurrentAtomicDomain;
class Archive;

// this is the map used internally by the atomic domain
#include "../data_structures/MutableMap.h"
typedef MutableMap<uint64_t, ConcurrentAtom*> ConcurrentAtomMapType;

struct ConcurrentAtomNeighborhood
{
    ConcurrentAtomNeighborhood();
    ConcurrentAtomNeighborhood(ConcurrentAtom *l, ConcurrentAtom *c, ConcurrentAtom *r);
    bool hasLeft() const;
    bool hasRight() const;
    ConcurrentAtom *center;
    ConcurrentAtom *left;
    ConcurrentAtom *right;
};

struct ConcurrentAtom
{
public:
    ConcurrentAtom(uint64_t p, float m);
    uint64_t pos() const;
    float mass() const;
    void updateMass(float newMass);
    friend Archive& operator<<(Archive& ar, const ConcurrentAtom &a);
    friend Archive& operator>>(Archive& ar, ConcurrentAtom &a);
//private: // TODO
    // only the atomic domain can change the position of an atom, since it is
    // responsible for keeping them ordered
    friend class ConcurrentAtomicDomain;
    void updatePos(uint64_t newPos);
    void setLeft(ConcurrentAtom *atom);
    void setRight(ConcurrentAtom *atom);
    void setIndex(unsigned index);
    void setIterator(ConcurrentAtomMapType::iterator it);    
    bool hasLeft() const;
    bool hasRight() const;
    ConcurrentAtom* left() const;
    ConcurrentAtom* right() const;
    unsigned index() const;
    ConcurrentAtomMapType::iterator iterator() const;

    uint64_t mPos;
    ConcurrentAtom *mLeft;
    ConcurrentAtom *mRight;
    ConcurrentAtomMapType::iterator mIterator; // iterator to position in map
    unsigned mIndex; // storing the index allows vector lookup once found in map
    float mMass;
};

#endif // __COGAPS_CONCURRENT_ATOM_H__