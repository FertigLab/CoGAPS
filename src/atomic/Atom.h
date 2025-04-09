#ifndef __COGAPS_ATOM_H__
#define __COGAPS_ATOM_H__

#include <map>

struct Atom;
class Archive;
class AtomicDomain;

// this is the map used internally by the atomic domain
// #include "../data_structures/MutableMap.h" -- it is unsafe
typedef std::map<uint64_t, size_t> AtomMapType;

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
    // only the atomic domain can change the position of an atom, since it is
    // responsible for keeping them ordered
    friend class AtomicDomain;
    bool hasLeft() const;
    bool hasRight() const;
    size_t leftIndex() const;
    size_t rightIndex() const;
    size_t index() const;
    AtomMapType::iterator iterator() const;

private:
    void updatePos(uint64_t newPos);
    void setLeftIndex(size_t index);
    void setRightIndex(size_t index);
    void unsetLeftIndex();
    void unsetRightIndex();
    void setIndex(size_t index);
    void setIterator(AtomMapType::iterator it);    

    AtomMapType::iterator mIterator; // iterator to position in map
    uint64_t mPos; // position of the atom
    bool mHasRight,mHasLeft; //are the following two defined
    std::size_t mLeftIndex; // index of left neighbor in the atomic storage
    std::size_t mRightIndex; // index of right neighbor in the atomic storage 
    std::size_t mIndex; // storing the index allows vector lookup once found in map
    float mMass; // mass of the atom
};

#endif // __COGAPS_ATOM_H__
