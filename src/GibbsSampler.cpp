#include "GibbsSampler.h"
#include "math/Algorithms.h"
#include "math/SIMD.h"

static float getDeltaLL(AlphaParameters alpha, float mass)
{
    return mass * (alpha.s_mu - alpha.s * mass / 2.f);
}

static OptionalFloat gibbsMass(AlphaParameters alpha, float a, float b,
float lambda, GapsRng *rng)
{
    if (alpha.s > gaps::epsilon)
    {
        float mean = (alpha.s_mu - lambda) / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        return rng->truncNormal(a, b, mean, sd);
    }
    return OptionalFloat();
}

unsigned GibbsSampler::dataRows() const
{
    return mDMatrix.nRow();
}

unsigned GibbsSampler::dataCols() const
{
    return mDMatrix.nCol();
}

void GibbsSampler::setSparsity(float alpha, float maxGibbsMass, bool singleCell)
{
    float meanD = singleCell ? gaps::algo::nonZeroMean(mDMatrix) :
        gaps::algo::mean(mDMatrix);

    mAlpha = alpha;
    mLambda = alpha * std::sqrt(mNumPatterns / meanD);
    mQueue.setAlpha(alpha);
    mQueue.setLambda(mLambda);
    mMaxGibbsMass = maxGibbsMass / mLambda;
}

void GibbsSampler::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

void GibbsSampler::setMatrix(const Matrix &mat)
{   
    GAPS_ASSERT(mat.nRow() == mMatrix.nRow());
    GAPS_ASSERT(mat.nCol() == mMatrix.nCol());

    mMatrix = mat;
}

float GibbsSampler::chi2() const
{   
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}
  
uint64_t GibbsSampler::nAtoms() const
{   
    return mDomain.size();
}

void GibbsSampler::recalculateAPMatrix()
{
    mAPMatrix = gaps::algo::matrixMultiplication(*mOtherMatrix, mMatrix);
}

void GibbsSampler::sync(const GibbsSampler &sampler, unsigned nThreads)
{
    mOtherMatrix = &(sampler.mMatrix);
    gaps::algo::copyTranspose(&mAPMatrix, sampler.mAPMatrix, nThreads);
}

void GibbsSampler::update(unsigned nSteps, unsigned nCores)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        // populate queue, prepare domain for this queue
        mQueue.populate(mDomain, nSteps - n);
        n += mQueue.size();
        
        // process all proposed updates
        #pragma omp parallel for num_threads(nCores)
        for (unsigned i = 0; i < mQueue.size(); ++i)
        {
            processProposal(mQueue[i]);
        }
        mQueue.clear();
    }

    GAPS_ASSERT(internallyConsistent());
    GAPS_ASSERT(mDomain.isSorted());
}

void GibbsSampler::processProposal(const AtomicProposal &prop)
{
    switch (prop.type)
    {
        case 'B':
            birth(prop);
            break;
        case 'D':
            death(prop);
            break;
        case 'M':
            move(prop);
            break;
        case 'E':
            exchange(prop);
            break;
    }
}

// add an atom at a random position, calculate mass either with an
// exponential distribution or with the gibbs mass distribution
void GibbsSampler::birth(const AtomicProposal &prop)
{
    // try to get mass using gibbs, resort to exponential if needed
    OptionalFloat mass = canUseGibbs(prop.c1)
        ? gibbsMass(alphaParameters(prop.r1, prop.c1) * mAnnealingTemp, 0.f,
            mMaxGibbsMass, mLambda, &(prop.rng))
        : prop.rng.exponential(mLambda);

    // accept mass as long as gibbs succeded or it's non-zero
    if (mass.hasValue() && mass.value() >= gaps::epsilon)
    {
        mQueue.acceptBirth();
        prop.atom1->mass = mass.value();
        changeMatrix(prop.r1, prop.c1, mass.value());
    }
    else
    {
        mQueue.rejectBirth();
        mDomain.erase(prop.atom1->pos);
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass distribution
void GibbsSampler::death(const AtomicProposal &prop)
{
    // calculate alpha parameters assuming atom dies
    AlphaParameters alpha = alphaParametersWithChange(prop.r1, prop.c1,
        -prop.atom1->mass);

    float rebirthMass = prop.atom1->mass;
    if (canUseGibbs(prop.c1))
    {
        OptionalFloat gMass = gibbsMass(alpha * mAnnealingTemp, 0.f,
            mMaxGibbsMass, mLambda, &(prop.rng));
        if (gMass.hasValue())
        {
            rebirthMass = gMass.value();
        }
    }

    // accept/reject rebirth
    float deltaLL = getDeltaLL(alpha, rebirthMass) * mAnnealingTemp;
    if (std::log(prop.rng.uniform()) < deltaLL) // accept rebirth
    {
        mQueue.rejectDeath();
        if (rebirthMass != prop.atom1->mass)
        {
            safelyChangeMatrix(prop.r1, prop.c1, rebirthMass - prop.atom1->mass);
        }
        prop.atom1->mass = rebirthMass;
    }
    else // reject rebirth
    {
        mQueue.acceptDeath();
        safelyChangeMatrix(prop.r1, prop.c1, -prop.atom1->mass);
        mDomain.erase(prop.atom1->pos);
    }
}

// move mass from src to dest in the atomic domain
void GibbsSampler::move(const AtomicProposal &prop)
{
    AlphaParameters alpha = alphaParameters(prop.r1, prop.c1, prop.r2, prop.c2);
    if (std::log(prop.rng.uniform()) < getDeltaLL(alpha, -prop.atom1->mass) * mAnnealingTemp)
    {
        prop.atom1->pos = prop.pos;
        safelyChangeMatrix(prop.r1, prop.c1, -prop.atom1->mass);
        changeMatrix(prop.r2, prop.c2, prop.atom1->mass);
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small
void GibbsSampler::exchange(const AtomicProposal &prop)
{
    // attempt gibbs distribution exchange
    AlphaParameters alpha = alphaParameters(prop.r1, prop.c1, prop.r2, prop.c2);
    if (canUseGibbs(prop.c1, prop.c2))
    {
        OptionalFloat gMass = gibbsMass(alpha * mAnnealingTemp,
            -prop.atom1->mass, prop.atom2->mass, 0.f, &(prop.rng));
        if (gMass.hasValue())
        {
            acceptExchange(prop, gMass.value());
            return;
        }
    }

    // resort to metropolis-hastings if gibbs fails
    exchangeUsingMetropolisHastings(prop, alpha);
}

void GibbsSampler::exchangeUsingMetropolisHastings(const AtomicProposal &prop,
AlphaParameters alpha)
{
    // compute amount of mass to be exchanged
    float totalMass = prop.atom1->mass + prop.atom2->mass;
    float newMass = prop.rng.truncGammaUpper(totalMass, 1.f / mLambda);

    // compute amount to change atom1 by - always change larger mass to newMass
    float delta = (prop.atom1->mass > prop.atom2->mass)
        ? newMass - prop.atom1->mass
        : prop.atom2->mass - newMass;

    // choose mass for priorLL calculation
    float oldMass = (2.f * newMass > totalMass)
        ? gaps::max(prop.atom1->mass, prop.atom2->mass)
        : gaps::min(prop.atom1->mass, prop.atom2->mass);

    // calculate priorLL
    float pNew = gaps::d_gamma(newMass, 2.f, 1.f / mLambda);
    float pOld = gaps::d_gamma(oldMass, 2.f, 1.f / mLambda);
    float priorLL = (pNew == 0.f) ? 1.f : pOld / pNew;

    // accept/reject
    float deltaLL = getDeltaLL(alpha, delta) * mAnnealingTemp;
    if (priorLL == 0.f || std::log(prop.rng.uniform() * priorLL) < deltaLL)
    {
        acceptExchange(prop, delta);
        return;
    }
}

// helper function for exchange step
void GibbsSampler::acceptExchange(const AtomicProposal &prop, float delta)
{
    if (prop.atom1->mass + delta > gaps::epsilon && prop.atom2->mass - delta > gaps::epsilon)
    {
        prop.atom1->mass += delta;
        prop.atom2->mass -= delta;

        changeMatrix(prop.r1, prop.c1, delta);
        changeMatrix(prop.r2, prop.c2, -delta);
    }
}

// here mass + delta is guaranteed to be positive
void GibbsSampler::changeMatrix(unsigned row, unsigned col, float delta)
{
    mMatrix(row, col) += delta;
    updateAPMatrix(row, col, delta);

    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
}

// delta could be negative, this is needed to prevent negative values in matrix
void GibbsSampler::safelyChangeMatrix(unsigned row, unsigned col, float delta)
{
    float newVal = gaps::max(mMatrix(row, col) + delta, 0.f);
    updateAPMatrix(row, col, newVal - mMatrix(row, col));
    mMatrix(row, col) = newVal;
}

void GibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->colPtr(col);
    float *ap = mAPMatrix.colPtr(row);
    unsigned size = mAPMatrix.nRow();

    gaps::simd::PackedFloat pOther, pAP;
    gaps::simd::Index i(0);
    gaps::simd::PackedFloat pDelta(delta);
    for (; i <= size - i.increment(); ++i)
    {
        pOther.load(other + i);
        pAP.load(ap + i);
        pAP += pDelta * pOther;
        pAP.store(ap + i);
    }

    for (unsigned j = i.value(); j < size; ++j)
    {
        ap[j] += delta * other[j];
    }
}

bool GibbsSampler::canUseGibbs(unsigned col) const
{
    return !gaps::algo::isVectorZero(mOtherMatrix->colPtr(col),
        mOtherMatrix->nRow());
}

bool GibbsSampler::canUseGibbs(unsigned c1, unsigned c2) const
{
    return canUseGibbs(c1) || canUseGibbs(c2);
}

AlphaParameters GibbsSampler::alphaParameters(unsigned row, unsigned col)
{
    return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(row),
        mSMatrix.colPtr(row), mAPMatrix.colPtr(row), mOtherMatrix->colPtr(col));
}

AlphaParameters GibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(r1),
            mSMatrix.colPtr(r1), mAPMatrix.colPtr(r1), mOtherMatrix->colPtr(c1),
            mOtherMatrix->colPtr(c2));
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

AlphaParameters GibbsSampler::alphaParametersWithChange(unsigned row,
unsigned col, float ch)
{
    return gaps::algo::alphaParametersWithChange(mDMatrix.nRow(),
        mDMatrix.colPtr(row), mSMatrix.colPtr(row), mAPMatrix.colPtr(row),
        mOtherMatrix->colPtr(col), ch);
}

Archive& operator<<(Archive &ar, GibbsSampler &s)
{
    ar << s.mMatrix << s.mDomain << s.mQueue << s.mAlpha
        << s.mLambda << s.mMaxGibbsMass << s.mAnnealingTemp << s.mNumPatterns
        << s.mNumBins << s.mBinLength << s.mDomainLength;
    return ar;
}

Archive& operator>>(Archive &ar, GibbsSampler &s)
{
    ar >> s.mMatrix >> s.mDomain >> s.mQueue >> s.mAlpha
        >> s.mLambda >> s.mMaxGibbsMass >> s.mAnnealingTemp >> s.mNumPatterns
        >> s.mNumBins >> s.mBinLength >> s.mDomainLength;
    return ar;
}

#ifdef GAPS_DEBUG
bool GibbsSampler::internallyConsistent()
{
    if (mDomain.begin() == mDomain.end())
    {
        float sum = gaps::algo::sum(mMatrix);
        if (sum != 0.f)
        {
            gaps_printf("non-zero matrix (%f) with zero domain\n", sum);
        }        
        return sum == 0.f;
    }

    std::vector<Atom*>::iterator it = mDomain.begin();
    float current = (*it)->mass;
    unsigned row = ((*it)->pos / mBinLength) / mNumPatterns;
    unsigned col = ((*it)->pos / mBinLength) % mNumPatterns;
    ++it;
    
    for (; it != mDomain.end(); ++it)
    {
        if (((*it)->pos / mBinLength) / mNumPatterns != row
        || ((*it)->pos / mBinLength) % mNumPatterns != col)
        {
            if (std::abs(current - mMatrix(row, col)) > 0.1f)
            {
                Rprintf("mass difference detected at row %lu, column %lu: %f %f\n",
                    row, col, current, mMatrix(row, col));
                return false;
            }
            current = (*it)->mass;
            row = ((*it)->pos / mBinLength) / mNumPatterns;
            col = ((*it)->pos / mBinLength) % mNumPatterns;
        }
        else
        {
            current += (*it)->mass;
        }
    }   
    return true;  
}
#endif // GAPS_DEBUG
