#ifndef __COGAPS_ATOM_H__
#define __COGAPS_ATOM_H__

#include "../utils/Archive.h"

struct Atom;
class AtomicDomain;
class ConcurrentAtomicDomain;

// this is the map used internally by the atomic domain
#include "../data_structures/MutableMap.h"
typedef MutableMap<uint64_t, Atom*> AtomMapType;

struct Atom
{
public:

    Atom(uint64_t p, float m);

    uint64_t pos() const;
    float mass() const;

    bool hasLeft() const;
    bool hasRight() const;
    Atom* left() const;
    Atom* right() const;

    void updateMass(float newMass);

    friend Archive& operator<<(Archive& ar, const Atom &a);
    friend Archive& operator>>(Archive& ar, Atom &a);

private:

    // only the atomic domain can change the position of an atom, since it is
    // responsible for keeping them ordered
    friend class AtomicDomain;
    friend class ConcurrentAtomicDomain;
    void updatePos(uint64_t newPos);

    void setLeft(Atom* atom);
    void setRight(Atom* atom);
    void setIndex(unsigned index);
    void setIterator(AtomMapType::iterator it);    

    unsigned index() const;
    AtomMapType::iterator iterator() const;

    Atom* mLeft;
    Atom* mRight;

    uint64_t mPos;

    AtomMapType::iterator mIterator; // iterator to position in map

    unsigned mIndex; // storing the index allows vector lookup once found in map
    float mMass;
};

#endif // __COGAPS_ATOM_H__