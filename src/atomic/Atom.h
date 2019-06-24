#ifndef __COGAPS_ATOM_H__
#define __COGAPS_ATOM_H__

struct Atom;
class Archive;
class AtomicDomain;

// this is the map used internally by the atomic domain
#include "../data_structures/MutableMap.h"
typedef MutableMap<uint64_t, unsigned> AtomMapType;

struct AtomNeighborhood
{
    AtomNeighborhood();
    AtomNeighborhood(Atom *l, Atom *c, Atom *r);
    bool hasLeft() const;
    bool hasRight() const;

    Atom* center;
    Atom* left;
    Atom* right;
};

struct Atom
{
public:
    Atom(uint64_t p, float m);
    uint64_t pos() const;
    float mass() const;
    void updateMass(float newMass);
    friend Archive& operator<<(Archive& ar, const Atom &a);
    friend Archive& operator>>(Archive& ar, Atom &a);
private:
    // only the atomic domain can change the position of an atom, since it is
    // responsible for keeping them ordered
    friend class AtomicDomain;
    void updatePos(uint64_t newPos);
    void setLeftIndex(int index);
    void setRightIndex(int index);
    void setIndex(int index);
    void setIterator(AtomMapType::iterator it);    
    bool hasLeft() const;
    bool hasRight() const;
    int leftIndex() const;
    int rightIndex() const;
    int index() const;
    AtomMapType::iterator iterator() const;

    AtomMapType::iterator mIterator; // iterator to position in map
    uint64_t mPos; // position of the atom
    int mLeftIndex; // index of left neighbor
    int mRightIndex; // index of right neighbor
    int mIndex; // storing the index allows vector lookup once found in map
    float mMass; // mass of the atom
};

#endif // __COGAPS_ATOM_H__