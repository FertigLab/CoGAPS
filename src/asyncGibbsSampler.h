// CRTP
template <class T, class M, class N>
class GibbsSampler
{
private:

    N* mDMatrix;
    N* mSMatrix;
    M* mMatrix;
    N* mAPMatrix; 

    AtomicQueue mQueue;
    AtomicDomain mDomain;
    float mLambda;

    unsigned mNumRows;
    unsigned mNumCols;
    unsigned mBinSize;

public:

    T* impl()
    {
        return static_cast<T*>(this);
    }

    void update(unsigned nSteps)
    {
        unsigned i = 0;
        while (i < nSteps)
        {
            // want this to be as quick as possible - otherwise there would be
            // a large speed up to making this run concurrently along with the
            // processProposal jobs, but that is much, much more complicated
            // to implement
            mQueue.populate(nSteps - mQueue.size() - i); 

            unsigned nJobs = std::min(mQueue.size(), nCores);
            for (unsigned i = 0; i < nJobs; ++i) // can be run in parallel
            {
                processProposal(mQueue[i]);
            }
            mQueue.clear(nJobs);
            i += nJobs;
        }
    }

    void processProposal(const AtomicProposal &prop)
    {
        switch (prop.type)
        {
            case 'D': death(prop.pos1, prop.mass1); break;
            case 'B': birth(prop.pos1, prop.mass1); break;
            case 'M': move(prop.pos1, prop.mass1, prop.pos2); break;
            case 'E': exchange(prop.pos1, prop.mass1, prop.pos2, prop.mass2); break;
        }
    }

    void birth(uint64_t pos, float mass)
    {
        unsigned row = impl->getRow(pos);
        unsigned col = impl->getCol(pos);
        float mass = impl()->canUseGibbs(row, col) ? gibbsMass(row, col)
            : gaps::random::exponential(mLambda);
        mDomain.addAtom(pos, mass);
        mMatrix(row, col) += mass;
        impl()->updateAPMatrix(row, col, mass);
    }

    void death(uint64_t pos, unsigned row, unsigned col)
    {
        float mass = mDomain.at(pos);
        mMatrix(row, col) += -mass;
        impl()->updateAPMatrix(row, col, -mass);

        float newMass = impl()->canUseGibbs(row, col) ? gibbsMass(row, col)
            : mass;
        float deltaLL = impl()->computeDeltaLL(row, col, newMass);

        if (deltaLL * mAnnealingTemp >= std::log(gaps::random::uniform()))
        {
            newMass = mDomain.updateAtomMass(pos, newMass - mass);
            mMatrix(row, col) += newMass;
            impl()->updateAPMatrix(row, col, newMass);
            if (mDomain.at(pos))
            {
                ++mQueue.minAtoms;
            }
        }
        else
        {
            mDomain.deleteAtom(pos);
            mQueue.maxAtoms--;
        }
    }

    void move(uint64_t p1, unsigned r1, unsigned c1, uint64_t p2, unsigned r2, unsigned c2)
    {
        float mass = mDomiain.at(p1);
        if (r1 == r2 && c1 == c2)
        {
            mDomain.deleteAtom(p1);
            mDomain.addAtom(p2, mass);
        }
        else
        {
            if (deltaLL * mAnnealingTemp >= std::log(gaps::random::uniform()))
            {
                mDomain.deleteAtom(p1);
                mDomain.addAtom(p2, mass);
                mMatrix(r1, c1) += -mass;
                mMatrix(r2, c2) += mass;
                impl()->updateAPMatrix(r1, c1, -mass);
                impl()->updateAPMatrix(r2, c2, mass);
            }
        }

    }

    void exchange(uint64_t p1, unsigned r1, unsigned c1, uint64_t p2, unsigned r2, unsigned c2)
    {
        float mass1 = mDomiain.at(p1);
        float mass2 = mDomiain.at(p2);
        float newMass = gaps::random::inverseGammaSample(0.f, mass1 + mass2, 2.f, 1.f / mLambda);
        float delta1 = mass1 > mass2 ? newMass - mass1 : mass2 - newMass;
        float delta2 = mass2 > mass1 ? newMass - mass2 : mass1 - newMass;
        if (r1 == r2 && c1 == c2)
        {
            mDomain.updateAtomMass(p1, delta1);
            mDomain.updateAtomMass(p2, delta2);
        }
        else
        {
            // all the exchange code
        }
    }

    float gibbsMass(unsigned row, unsigned col)
    {        
        AlphaParameters alpha = impl()->alphaParameters(row, col);
        alpha.s *= mAnnealingTemp / 2.f;
        alpha.su *= mAnnealingTemp / 2.f;
        float mean  = (2.f * alpha.su - mLambda) / (2.f * alpha.s);
        float sd = 1.f / std::sqrt(2.f * alpha.s);

        float plower = gaps::random::p_norm(0.f, mean, sd);
        
        float newMass = death ? mass : 0.f;
        if (plower < 1.f && alpha.s > 0.00001f)
        {
            newMass = gaps::random::inverseNormSample(plower, 1.f, mean, sd);
        }
        return std::max(0.f, std::min(newMax, mMaxGibbsMass));
    }
};

class AmplitudeGibbsSampler : public GibbsSampler<AmplitudeGibbsSampler, ColMatrix, RowMatrix>
{
public:

    bool canUseGibbs(unsigned row, unsigned col)
    {
        return !gaps::algo::isRowZero(mOtherMatrix, col);
    }

    bool canUseGibbs(unsigned row1, unsigned col1, unsigned row2, unsigned col2)
    {
        return !gaps::algo::isRowZero(mOtherMatrix, col1)
            && !gaps::algo::isRowZero(mOtherMatrix, col2);
    }

    void updateAPMatrix(unsigned row, unsigned col, float delta)
    {
        for (unsigned j = 0; j < mAPMatrix.nCol(); ++j)
        {
            mAPMatrix(row,j) += delta * mOtherMatrix(j,col);
        }
    }

    unsigned getRow(uint64_t pos) const
    {
        return std::min(pos / (mBinSize * mNumCols), mNumRows);
    }
    
    unsigned getCol(uint64_t pos) const
    {
        return std::min((pos / mBinSize) % mNumCols, mNumCols);
    }
};

class PatternGibbsSampler : public GibbsSampler<PatternGibbsSampler, RowMatrix, ColMatrix>
{
public:

    bool canUseGibbs(unsigned row, unsigned col)
    {
        return !gaps::algo::isColZero(mOtherMatrix, row);
    }

    bool canUseGibbs(unsigned row1, unsigned col1, unsigned row2, unsigned col2)
    {
        return !gaps::algo::isColZero(mOtherMatrix, row1)
            && !gaps::algo::isColZero(mOtherMatrix, row2);
    }

    void updateAPMatrix(unsigned row, unsigned col, float delta)
    {
        for (unsigned i = 0; i < mAPMatrix.nRow(); ++i)
        {
            mAPMatrix(i,col) += delta * mOtherMatrix(i,row);
        }
    }

    unsigned getRow(uint64_t pos) const
    {
        return std::min((pos / mBinSize) % mNumRows, mNumRows);
    }
    
    unsigned getCol(uint64_t pos) const
    {
        return std::min(pos / (mBinSize * mNumRows), mNumCols);
    }
};
