#include "GibbsSampler.h"
#include "math/SIMD.h"

static float getDeltaLL(AlphaParameters alpha, float mass)
{
    return mass * (alpha.su - alpha.s * mass / 2.f);
}

unsigned GibbsSampler::dataRows() const
{
    return mDMatrix.nRow();
}

unsigned GibbsSampler::dataCols() const
{
    return mDMatrix.nCol();
}

void GibbsSampler::setSparsity(float alpha, bool singleCell)
{
    float meanD = singleCell ? gaps::algo::nonZeroMean(mDMatrix) :
        gaps::algo::mean(mDMatrix);

    mAlpha = alpha;
    mLambda = alpha * std::sqrt(mNumPatterns / meanD);
}

void GibbsSampler::setMaxGibbsMass(float max)
{
    mMaxGibbsMass = max;
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

void GibbsSampler::sync(const GibbsSampler &sampler)
{
    mOtherMatrix = &(sampler.mMatrix);
    gaps::algo::copyTranspose(&mAPMatrix, sampler.mAPMatrix);
}

void GibbsSampler::update(unsigned nSteps, unsigned nCores)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        makeAndProcessProposal();
        ++n;
    }
}

void GibbsSampler::makeAndProcessProposal()
{
    // always birth when no atoms exist
    if (mDomain.size() == 0)
    {
        return birth();
    }

    float bdProb = mDomain.size() < 2 ? 0.6667f : 0.5f;

    float u1 = mPropRng.uniform();
    float u2 = mPropRng.uniform();

    float lowerBound = deathProb(mDomain.size());

    if (u1 <= bdProb)
    {
        return u2 < lowerBound ? death() : birth();
    }
    return (u1 < 0.75f || mDomain.size() < 2) ? move() : exchange();
}

float GibbsSampler::deathProb(uint64_t nAtoms) const
{
    double size = static_cast<double>(mDomainLength);
    double term1 = (size - static_cast<double>(nAtoms)) / size;
    double term2 = mAlpha * static_cast<double>(mNumBins) * term1;
    return static_cast<double>(nAtoms) / (static_cast<double>(nAtoms) + term2);
}

// add an atom at a random position, calculate mass either with an
// exponential distribution or with the gibbs mass distribution
void GibbsSampler::birth()
{
    uint64_t pos = mDomain.randomFreePosition();
    unsigned row = getRow(pos);
    unsigned col = getCol(pos);
    GapsRng rng;

    // calculate proposed mass
    float mass = canUseGibbs(col)
        ? gibbsMass(alphaParameters(row, col), &rng).value
        : rng.exponential(mLambda);

    // accept mass as long as it's non-zero
    if (mass >= gaps::epsilon)
    {
        mDomain.insert(pos, mass);
        mMatrix(row, col) += mass;
        updateAPMatrix(row, col, mass);
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass calculation
void GibbsSampler::death()
{
    // get random atom
    Atom* atom = mDomain.randomAtom();
    unsigned row = getRow(atom->pos);
    unsigned col = getCol(atom->pos);
    GapsRng rng;

    // kill off atom
    safelyChangeMatrix(row, col, -1.f * atom->mass);

    // calculate rebirth mass
    AlphaParameters alpha = alphaParameters(row, col);
    float rebirthMass = atom->mass;
    if (canUseGibbs(col))
    {
        OptionalFloat gMass = gibbsMass(alpha, &rng);
        if (gMass.hasValue)
        {
            rebirthMass = gMass.value;
        }
    }

    // accept/reject rebirth
    if (getDeltaLL(alpha, rebirthMass) * mAnnealingTemp >= std::log(rng.uniform()))
    {
        atom->mass = rebirthMass;
        mMatrix(row, col) += rebirthMass;
        updateAPMatrix(row, col, rebirthMass);
    }
    else
    {
        mDomain.erase(atom->pos);
    }
}

// move mass from src to dest in the atomic domain
void GibbsSampler::move()
{
    AtomNeighborhood hood = mDomain.randomAtomWithNeighbors();
    uint64_t lbound = hood.hasLeft() ? hood.left->pos : 0;
    uint64_t rbound = hood.hasRight() ? hood.right->pos : mDomainLength;
    uint64_t newLocation = mPropRng.uniform64(lbound + 1, rbound - 1);

    unsigned r1 = getRow(hood.center->pos);
    unsigned c1 = getCol(hood.center->pos);
    unsigned r2 = getRow(newLocation);
    unsigned c2 = getCol(newLocation);

    if (r1 != r2 || c1 != c2)
    {
        GapsRng rng;
        AlphaParameters alpha = alphaParameters(r1, c1, r2, c2);
        float deltaLL = getDeltaLL(alpha, -1.f * hood.center->mass);
        if (deltaLL * mAnnealingTemp > std::log(rng.uniform()))
        {
            hood.center->pos = newLocation;
            safelyChangeMatrix(r1, c1, -1.f * hood.center->mass);
            mMatrix(r2, c2) += hood.center->mass;
            updateAPMatrix(r2, c2, hood.center->mass);
        }
    }
    else
    {
        hood.center->pos = newLocation;
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small after
// the exchange
void GibbsSampler::exchange()
{
    AtomNeighborhood hood = mDomain.randomAtomWithRightNeighbor();
    Atom* a1 = hood.center;
    Atom* a2 = hood.hasRight() ? hood.right : mDomain.front();

    unsigned r1 = getRow(a1->pos);
    unsigned c1 = getCol(a1->pos);
    unsigned r2 = getRow(a2->pos);
    unsigned c2 = getCol(a2->pos);

    if (r1 != r2 || c1 != c2)
    {
        GapsRng rng;

        float m1 = a1->mass;
        float m2 = a2->mass;

        AlphaParameters alpha = alphaParameters(r1, c1, r2, c2);
        if (canUseGibbs(c1, c2))
        {
            OptionalFloat gMass = gibbsMass(alpha, m1, m2, &rng);
            if (gMass.hasValue)
            {
                acceptExchange(a1, a2, gMass.value, r1, c1, r2, c2);
                return;
            }
        }

        float newMass = rng.truncGammaUpper(m1 + m2, 2.f, 1.f / mLambda);

        float delta = m1 > m2 ? newMass - m1 : m2 - newMass; // change larger mass
        float pOldMass = 2.f * newMass > m1 + m2 ? gaps::max(m1, m2) : gaps::min(m1, m2);

        float pNew = gaps::d_gamma(newMass, 2.f, 1.f / mLambda);
        float pOld = gaps::d_gamma(pOldMass, 2.f, 1.f / mLambda);

        float deltaLL = getDeltaLL(alpha, delta);
        float priorLL = (pNew == 0.f) ? 1.f : pOld / pNew;
        if (priorLL == 0.f || std::log(rng.uniform() * priorLL) < deltaLL * mAnnealingTemp)
        {
            acceptExchange(a1, a2, delta, r1, c1, r2, c2);
            return;
        }
    }
}

// helper function for exchange step
void GibbsSampler::acceptExchange(Atom *a1, Atom *a2, float delta,
unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    float d1 = updateAtomMass(a1, delta) ? delta : -a1->mass;
    float d2 = updateAtomMass(a2, -delta) ? -delta : -a2->mass;

    safelyChangeMatrix(r1, c1, d1);
    safelyChangeMatrix(r2, c2, d2);
}

// helper function for acceptExchange, used to conditionally update the mass
// at a single position, deleting the atom if the new mass is too small
bool GibbsSampler::updateAtomMass(Atom *atom, float delta)
{
    if (atom->mass + delta < gaps::epsilon)
    {
        mDomain.erase(atom->pos);
        return false;
    }
    atom->mass += delta;
    return true;
}

OptionalFloat GibbsSampler::gibbsMass(AlphaParameters alpha,
GapsRng *rng)
{        
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = (alpha.su - mLambda) / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::p_norm(0.f, mean, sd);

        if (pLower < 1.f)
        {
            float m = rng->inverseNormSample(pLower, 1.f, mean, sd);
            float gMass = gaps::min(m, mMaxGibbsMass / mLambda);
            if (gMass >= gaps::epsilon)
            {
                return OptionalFloat(gMass);
            }
        }
    }
    return OptionalFloat();
}

OptionalFloat GibbsSampler::gibbsMass(AlphaParameters alpha,
float m1, float m2, GapsRng *rng)
{
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = alpha.su / alpha.s; // lambda cancels out
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::p_norm(-m1, mean, sd);
        float pUpper = gaps::p_norm(m2, mean, sd);

        if (!(pLower >  0.95f || pUpper < 0.05f))
        {
            float delta = rng->inverseNormSample(pLower, pUpper, mean, sd);
            float gMass = gaps::min(gaps::max(-m1, delta), m2); // conserve mass
            return OptionalFloat(gMass);
        }
    }
    return OptionalFloat();
}

// needed to prevent negative values in matrix
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

    gaps::simd::packedFloat pOther, pAP;
    gaps::simd::Index i(0);
    gaps::simd::packedFloat pDelta(delta);
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

unsigned GibbsSampler::getRow(uint64_t pos) const
{
    return pos / (mBinSize * mNumPatterns); // nCol == nPatterns
}

unsigned GibbsSampler::getCol(uint64_t pos) const
{
    return (pos / mBinSize) % mNumPatterns; // nCol == nPatterns
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

Archive& operator<<(Archive &ar, GibbsSampler &s)
{
    // TODO
    return ar;
}

Archive& operator>>(Archive &ar, GibbsSampler &s)
{
    // TODO
    return ar;
}

#ifdef GAPS_DEBUG
bool GibbsSampler::internallyConsistent()
{
    return true;
}
#endif // GAPS_DEBUG
