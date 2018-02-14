#ifndef __GAPS_ATOMIC_DOMAIN_H__
#define __GAPS_ATOMIC_DOMAIN_H__

// data structure that holds atoms
class AtomicDomain
{
private:

public:

    AtomicDomain();
    
    uint64_t randomAtomPosition();
    uint64_t randomFreePosition();

    float updateMass(uint64_t pos, float delta);

};

#endif


struct Atom
{
    uint64_t pos;
    float mass;

    Atom* left;
    Atom* right;
    
    Atom(uint64_t p, float m)
        : pos(p), mass(m), left(nullptr), right(nullptr)
    {}
};

void insertAtom(uint64_t p, float m)
{
    std::map<uint64_t, Atom>::const_iterator it, left, right;
    it = mAtoms.insert(std::pair<uint64_t, Atom>(p, Atom(p,m))).first;
    
    std::map<uint64_t, Atom>::const_iterator left(it), right(it);
    if (it != mAtoms.begin())
    {
        --left;
    }
    if (++it != mAtoms.end())
    {
        ++right;
    }
    it->left = &left;
    it->right = &right;
}



void removeAtom(const Atom &atom)
{
    atom.left.right = atom.right;
    atom.right.left = atom.left;
}
