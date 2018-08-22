#include "AtomicDomain.h"
#include "GapsAssert.h"

#include <vector>

// should only ever use this when copying a default constructed bucket
// in the vector of buckets used in the allocator
/*AtomBucket& AtomBucket::operator=(const AtomBucket& other)
{
    GAPS_ASSERT(other.mSize == 0);
    GAPS_ASSERT(other.mPrev == NULL);
    GAPS_ASSERT(other.mNext == NULL);
    GAPS_ASSERT(other.mOverflow == NULL);

    mSize = 0;
    mPrev = NULL;
    mNext = NULL;
    mOverflow = NULL;
}*/

/////////////////////// ALLOCATOR FOR OVERFLOW BUCKETS /////////////////////////

static const unsigned poolsize = 1024;

class AtomBucketPool
{
public:

    AtomBucketPool()
        : mPool(std::vector<AtomBucket>(poolsize, AtomBucket())), mCurrent(0)
    {}

    AtomBucket* create()
    {
        return &(mPool[mCurrent++]);
    }

    bool depleted()
    {
        return mCurrent == poolsize;
    }

private:

    std::vector<AtomBucket> mPool;
    unsigned mCurrent;
};

class AtomBucketAllocator
{
public:

    AtomBucketAllocator()
    {
        mAllPools.push_back(new AtomBucketPool());
    }

    ~AtomBucketAllocator()
    {
        for (unsigned i = 0; i < mAllPools.size(); ++i)
        {
            delete mAllPools[i];
        }
    }

    AtomBucket* create()
    {
        if (mAllPools.back()->depleted())
        {
            mAllPools.push_back(new AtomBucketPool());
        }
        return mAllPools.back()->create();
    }

private:

    std::vector<AtomBucketPool*> mAllPools;
};

static AtomBucket* createAtomBucket()
{
    //static AtomBucketAllocator allocator;
    //return allocator.create();
    return new AtomBucket();
}

////////////////////////////////// ATOM ////////////////////////////////////////

Atom::Atom() : pos(0), mass(0.f) {}

Atom::Atom(uint64_t p, float m) : pos(p), mass(m) {}

void Atom::operator=(Atom other)
{
    pos = other.pos;
    mass = other.mass;
}

Archive& operator<<(Archive &ar, Atom &a)
{
    ar << a.pos << a.mass;
    return ar;
}

Archive& operator>>(Archive &ar, Atom &a)
{
    ar >> a.pos >> a.mass;
    return ar;
}

//////////////////////////// ATOM NEIGHBORHOOD /////////////////////////////////


AtomNeighborhood::AtomNeighborhood(Atom *l, Atom *c, Atom *r)
    : center(c), left(l), right(r)
{}

bool AtomNeighborhood::hasLeft()
{
    return left != NULL;
}

bool AtomNeighborhood::hasRight()
{
    return right != NULL;
}

/////////////////////////////// ATOM BUCKET ////////////////////////////////////

// note much of this code is unrolled intentionally so that the most frequent
// cases are fast

AtomBucket::AtomBucket()
    : mSize(0), mOverflow(NULL), mPrev(NULL), mNext(NULL)
{}

unsigned AtomBucket::size() const
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());

    return mOverflow == NULL ? mSize : mSize + mOverflow->size();
}

bool AtomBucket::isEmpty() const
{
    return mSize == 0;
}

bool AtomBucket::contains(uint64_t pos) const
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(pos > 0);

    if (mOverflow != NULL && !mOverflow->isEmpty() && pos > mBucket[1].pos)
    {
        return mOverflow->contains(pos);
    }
    else
    {
        switch (mSize)
        {
            case 0:
                return false;
            case 1:
                return mBucket[0].pos == pos;
            case 2:
                return mBucket[0].pos == pos || mBucket[1].pos == pos;
        }
    }
}

Atom* AtomBucket::operator[](unsigned index)
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(index < size());

    if (index > 1) // must be overflowed
    {
        return mOverflow->operator[](index - 2);
    }
    return &(mBucket[index]);
}

void AtomBucket::insert(uint64_t pos, float mass)
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(pos > 0);

    if (mSize == 0)
    {
        mBucket[0] = Atom(pos, mass);
        ++mSize;
    }
    else if (mSize == 1)
    {
        if (pos < mBucket[0].pos)
        {
            mBucket[1] = mBucket[0];
            mBucket[0] = Atom(pos, mass);
        }
        else
        {
            mBucket[1] = Atom(pos, mass);
        }
        ++mSize;
    }
    else
    {
        // check if we need to allocate the overflow bucket
        if (mOverflow == NULL)
        {
            mOverflow = createAtomBucket();
            mOverflow->mPrev = this;
            mOverflow->mNext = mNext;
        }
        else if (mOverflow->isEmpty())
        {
            mOverflow->mPrev = this;
            mOverflow->mNext = mNext;
        }

        // push correct atom into overflow bucket
        if (pos > mBucket[1].pos)
        {
            return mOverflow->insert(pos, mass);
        }
        mOverflow->insert(mBucket[1].pos, mBucket[1].mass);

        // if inserting in this bucket, find correct position
        if (pos < mBucket[0].pos)
        {
            mBucket[1] = mBucket[0];
            mBucket[0] = Atom(pos, mass);
        }
        else
        {
            mBucket[1] = Atom(pos, mass);
        }
    }
}

// assumes pos is contained in this chain
void AtomBucket::erase(uint64_t pos)
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(pos > 0);
    GAPS_ASSERT(mSize > 0);
    GAPS_ASSERT(contains(pos));

    if (mSize == 1)
    {
        connectAdjacent();
        mSize = 0;
    }
    else if (mSize == 2)
    {
        // check if this position is in overflow bucket    
        if (pos > mBucket[1].pos)
        {
            return mOverflow->erase(pos);
        }

        // shift top position down if needed
        if (mBucket[0].pos == pos)
        {
            mBucket[0] = mBucket[1];
        }
        
        // pull first atom from overflow if it's there
        if (mOverflow != NULL && !mOverflow->isEmpty())
        {
            mBucket[1] = *(mOverflow->front());
            mOverflow->eraseFront();
        }
        else // just delete atom at last position
        {
            --mSize;
        }
    }
}

void AtomBucket::eraseFront()
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());

    GAPS_ASSERT(mSize > 0);

    if (mSize == 1)
    {
        connectAdjacent();
        mSize = 0;
    }
    else if (mSize == 2)
    {
        mBucket[0] = mBucket[1];
        
        // pull first atom from overflow if it's there
        if (mOverflow != NULL && !mOverflow->isEmpty())
        {
            mBucket[1] = *(mOverflow->front());
            mOverflow->eraseFront();
        }
        else
        {
            --mSize;
        }
    }
}

AtomNeighborhood AtomBucket::getNeighbors(unsigned index)
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(index < size());

    return AtomNeighborhood(getLeft(index), this->operator[](index), getRight(index));
}

AtomNeighborhood AtomBucket::getRightNeighbor(unsigned index)
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(index < size());

    return AtomNeighborhood(NULL, this->operator[](index), getRight(index));
}

// needs to propogate through overflow so that the last overflow will know
// where the next bucket is
void AtomBucket::setRightAdjacentBucket(AtomBucket *bucket)
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());

    mNext = bucket;
    if (mOverflow != NULL && !mOverflow->isEmpty())
    {
        mOverflow->setRightAdjacentBucket(bucket);
    }
}

void AtomBucket::setLeftAdjacentBucket(AtomBucket *bucket)
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(bucket == NULL || !bucket->isEmpty());

    mPrev = bucket;
}

void AtomBucket::connectAdjacent()
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());

    // be careful when connecting overflow buckets
    if (mNext != NULL && (mPrev == NULL || mPrev->mOverflow != this))
    {
        mNext->setLeftAdjacentBucket(mPrev);
    }
    if (mPrev != NULL)
    {
        mPrev->setRightAdjacentBucket(mNext);
    }

    mNext = NULL;
    mPrev = NULL;
}

Atom* AtomBucket::front()
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(mSize > 0);
    GAPS_ASSERT(mBucket[0].pos > 0);

    return &(mBucket[0]);
}

Atom* AtomBucket::back()
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(mSize > 0);

    if (mOverflow == NULL || (mOverflow != NULL && mOverflow->isEmpty()))
    {
        return &(mBucket[mSize - 1]);
    }
    return mOverflow->back();
}

Atom* AtomBucket::getLeft(unsigned index)
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(mSize > 0);
    GAPS_ASSERT(index < size());

    if (index == 0)
    {
        GAPS_ASSERT(mPrev == NULL || mPrev->mOverflow != this);
        return mPrev != NULL ? mPrev->back() : NULL; // mPrev can't be overflow here
    }
    else if (index < 3)
    {
        return &(mBucket[index - 1]);
    }
    else
    {
        GAPS_ASSERT(mOverflow != NULL);
        return mOverflow->getLeft(index - 2);
    }    
}

Atom* AtomBucket::getRight(unsigned index)
{
    GAPS_ASSERT(mPrev == NULL || !mPrev->isEmpty());
    GAPS_ASSERT(mNext == NULL || !mNext->isEmpty());
    GAPS_ASSERT(mSize > 0);
    GAPS_ASSERT(index < size());

    if (index == mSize - 1) // this is last atom in bucket
    {
        if (mOverflow != NULL && !mOverflow->isEmpty())
        {
            return mOverflow->front();
        }
        return mNext != NULL ? mNext->front() : NULL;
    }        
    else if (index == 0)
    {
        GAPS_ASSERT(mBucket[1].pos > 0);
        return &(mBucket[1]);
    }
    else
    {
        mOverflow->getRight(index - 2);
    }
}

Archive& operator<<(Archive &ar, AtomBucket &b)
{
    bool hasOverflow = (b.mOverflow != NULL);
    ar << b.mBucket[0] << b.mBucket[1] << b.mSize << hasOverflow;

    if (hasOverflow)
    {
        ar << *b.mOverflow;
    }
    return ar;
}

Archive& operator>>(Archive& ar, AtomBucket &b)
{
    bool hasOverflow = false;
    ar >> b.mBucket[0] >> b.mBucket[1] >> b.mSize >> hasOverflow;

    if (hasOverflow)
    {
        b.mOverflow = createAtomBucket();
        ar >> *b.mOverflow;
    }
    return ar;
}

////////////////////////////// ATOM HASH MAP ///////////////////////////////////

AtomHashMap::AtomHashMap(unsigned expectedNumAtoms)
    :
mExpectedNumAtoms(expectedNumAtoms), mLongestBucketSize(0), mSize(0),
mHashMap(std::vector<AtomBucket>(expectedNumAtoms, AtomBucket()))
{}

void AtomHashMap::setTotalLength(uint64_t length)
{
    GAPS_ASSERT(length % mExpectedNumAtoms == 0);
    mTotalLength = length;
    mBucketLength = mTotalLength / mExpectedNumAtoms;
}

// pos ranges from 1 to mTotalLength
unsigned AtomHashMap::hash(uint64_t pos) const
{
    return pos / mBucketLength;
}

Atom* AtomHashMap::front()
{
    return mHashMap[*(mFullBuckets.begin())].front();
}
    
HashMapIndex AtomHashMap::getRandomIndex() const
{
    GAPS_ASSERT(mSize > 0);
    while (true)
    {
        unsigned bucket = hash(mRng.uniform64(0, mTotalLength - 1));
        if (!mHashMap[bucket].isEmpty())
        {
            unsigned pos = mRng.uniform32(0, mLongestBucketSize - 1);
            if (pos < mHashMap[bucket].size())
            {
                return HashMapIndex(bucket, pos);
            }
        }
    }
}

Atom* AtomHashMap::randomAtom()
{
    HashMapIndex ndx = getRandomIndex();
    return mHashMap[ndx.bucket][ndx.index];
}

AtomNeighborhood AtomHashMap::randomAtomWithNeighbors()
{
    HashMapIndex ndx = getRandomIndex();
    return mHashMap[ndx.bucket].getNeighbors(ndx.index);
}

AtomNeighborhood AtomHashMap::randomAtomWithRightNeighbor()
{
    HashMapIndex ndx = getRandomIndex();
    return mHashMap[ndx.bucket].getRightNeighbor(ndx.index);
}

bool AtomHashMap::contains(uint64_t pos) const
{
    unsigned index = hash(pos);
    return !(mHashMap[index].isEmpty()) && mHashMap[index].contains(pos);
}

unsigned AtomHashMap::size() const
{
    return mSize;
}

// usually O(1), might be O(logN)
void AtomHashMap::insert(uint64_t pos, float mass)
{
    // hash position once
    unsigned index = hash(pos);

    // if inserting into an empty bucket, mark bucket as non-empty
    // and connect adjacent buckets O(logN)
    if (mHashMap[index].isEmpty())
    {
        mHashMap[index].insert(pos, mass);
        GAPS_ASSERT(mHashMap[index].contains(pos));
        std::set<unsigned>::iterator it = mFullBuckets.insert(index).first;
        if (it != mFullBuckets.begin())
        {
            --it;
            mHashMap[index].setLeftAdjacentBucket(&mHashMap[*it]);
            mHashMap[*it].setRightAdjacentBucket(&mHashMap[index]);
            ++it;
        }
        if (++it != mFullBuckets.end())
        {
            mHashMap[index].setRightAdjacentBucket(&mHashMap[*it]);
            mHashMap[*it].setLeftAdjacentBucket(&mHashMap[index]);
        }
    }
    else
    {
        mHashMap[index].insert(pos, mass);
        GAPS_ASSERT(mHashMap[index].contains(pos));
    }
    GAPS_ASSERT(mHashMap[index].contains(pos));

    // insert atom
    ++mSize;

    // check if this is now the longest bucket
    if (mHashMap[index].size() > mLongestBucketSize)
    {
        mWhichLongestBucket = index;
        ++mLongestBucketSize;
    }
}

// usually O(1), sometimes O(logN), rarely O(N)
void AtomHashMap::erase(uint64_t pos)
{
    // erase atom
    unsigned index = hash(pos);
    GAPS_ASSERT(!mHashMap[index].isEmpty());
    GAPS_ASSERT(mHashMap[index].contains(pos));
    --mSize;
    
    // mark bucket as empty if this was the last atom O(logN)
    mHashMap[index].erase(pos);
    if (mHashMap[index].isEmpty())
    {
        mFullBuckets.erase(index);
    }
    GAPS_ASSERT(!mHashMap[index].contains(pos));

    // if this atom was in the largest bucket, find the new largest bucket O(N)
    if (index == mWhichLongestBucket)
    {
        --mLongestBucketSize;
        std::set<unsigned>::iterator it = mFullBuckets.begin();
        for (; it != mFullBuckets.end(); ++it)
        {
            if (mHashMap[*it].size() > mLongestBucketSize)
            {
                mLongestBucketSize = mHashMap[*it].size();
                mWhichLongestBucket = *it;
            }
        }
    }
}

void AtomHashMap::move(uint64_t src, uint64_t dest, float mass)
{
    unsigned srcHash = hash(src);
    unsigned destHash = hash(dest);
    if (srcHash != destHash)
    {
        erase(src);
        insert(dest, mass);
    }
}

Archive& operator<<(Archive& ar, AtomHashMap &h)
{
    ar << h.mExpectedNumAtoms << h.mBucketLength << h.mTotalLength <<
        h.mWhichLongestBucket << h.mLongestBucketSize << h.mSize;

    ar << h.mFullBuckets.size();
    std::set<unsigned>::iterator it = h.mFullBuckets.begin();
    for (; it != h.mFullBuckets.end(); ++it)
    {
        ar << *it << h.mHashMap[*it];
    }
    return ar;
}

Archive& operator>>(Archive& ar, AtomHashMap &h)
{
    ar >> h.mExpectedNumAtoms >> h.mBucketLength >> h.mTotalLength >>
        h.mWhichLongestBucket >> h.mLongestBucketSize >> h.mSize;

    unsigned nBuckets = 0;
    ar >> nBuckets;

    for (unsigned i = 0; i < nBuckets; ++i)
    {
        unsigned thisBucket = 0;
        ar >> thisBucket;
        h.mFullBuckets.insert(thisBucket);
        ar >> h.mHashMap[thisBucket];
    }
    return ar;
}

////////////////////////////// ATOMIC DOMAIN ///////////////////////////////////

AtomicDomain::AtomicDomain(unsigned nBins) : mAtoms(nBins)
{
    uint64_t binLength = std::numeric_limits<uint64_t>::max() / nBins;
    mDomainLength = binLength * nBins;
    mAtoms.setTotalLength(mDomainLength);
}

void AtomicDomain::setDomainSize(uint64_t size)
{
    // do nothing
}

Atom* AtomicDomain::front()
{
    return mAtoms.front();
}

Atom* AtomicDomain::randomAtom()
{
    return mAtoms.randomAtom();
}

AtomNeighborhood AtomicDomain::randomAtomWithNeighbors()
{
    return mAtoms.randomAtomWithNeighbors();
}

AtomNeighborhood AtomicDomain::randomAtomWithRightNeighbor()
{
    return mAtoms.randomAtomWithRightNeighbor();
}

uint64_t AtomicDomain::randomFreePosition() const
{
    uint64_t pos = mRng.uniform64(1, mDomainLength);
    while (mAtoms.contains(pos))
    {
        pos = mRng.uniform64(1, mDomainLength);
    } 
    return pos;
}

uint64_t AtomicDomain::size() const
{
    return mAtoms.size();
}

void AtomicDomain::cacheInsert(uint64_t pos, float mass) const
{
    unsigned ndx = 0;
    #pragma omp critical(atomicInsert)
    {
        ndx = mInsertCacheIndex++;
    }
    mInsertCache[ndx] = Atom(pos, mass);
}

void AtomicDomain::cacheErase(uint64_t pos) const
{
    unsigned ndx = 0;
    #pragma omp critical(atomicErase)
    {
        ndx = mEraseCacheIndex++;
    }
    mEraseCache[ndx] = pos;
}

void AtomicDomain::cacheMove(uint64_t src, uint64_t dest, float mass) const
{
    unsigned ndx = 0;
    #pragma omp critical(atomicMove)
    {
        ndx = mMoveCacheIndex++;
    }
    mMoveCache[ndx] = MoveOperation(src, dest, mass);
}

void AtomicDomain::resetCache(unsigned n)
{
    mInsertCacheIndex = 0;
    mEraseCacheIndex = 0;
    mMoveCacheIndex = 0;
    mInsertCache.resize(n);
    mEraseCache.resize(n);
    mMoveCache.resize(n);
}

void AtomicDomain::flushCache()
{
    for (unsigned i = 0; i < mEraseCacheIndex; ++i)
    {
        mAtoms.erase(mEraseCache[i]);
    }

    for (unsigned i = 0; i < mInsertCacheIndex; ++i)
    {
        mAtoms.insert(mInsertCache[i].pos, mInsertCache[i].mass);
    }

    for (unsigned i = 0; i < mMoveCacheIndex; ++i)
    {
        mAtoms.move(mMoveCache[i].src, mMoveCache[i].dest, mMoveCache[i].mass);
    }

    mInsertCache.clear();
    mEraseCache.clear();
    mMoveCache.clear();
}

Archive& operator<<(Archive &ar, AtomicDomain &domain)
{
    ar << domain.mDomainLength << domain.mAtoms << domain.mRng;
}

Archive& operator>>(Archive &ar, AtomicDomain &domain)
{
    ar >> domain.mDomainLength >> domain.mAtoms << domain.mRng;
}
