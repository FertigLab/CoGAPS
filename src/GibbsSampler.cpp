#include "GibbsSampler.h"
#include "math/SIMD.h"

void GibbsSampler::sync(GibbsSampler &sampler)
{
    mOtherMatrix = &(sampler.mMatrix);
    gaps::algo::copyTranspose(&mAPMatrix, sampler.mAPMatrix);
}

void GibbsSampler::setSparsity(float alpha, bool singleCell)
{
    float meanD = singleCell ? gaps::algo::nonZeroMean(mDMatrix) :
        gaps::algo::mean(mDMatrix);

    unsigned nPatterns = mMatrix.nCol();

    mAlpha = alpha;
    mLambda = alpha * std::sqrt(nPatterns / meanD);
}

unsigned GibbsSampler::getRow(uint64_t pos) const
{
    return pos / (mBinSize * mNumCols);
}

unsigned GibbsSampler::getCol(uint64_t pos) const
{
    return (pos / mBinSize) % mNumCols;
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

float GibbsSampler::computeDeltaLL(unsigned row, unsigned col, float mass)
{
    AlphaParameters alpha = alphaParameters(row, col);
    return mass * (alpha.su - alpha.s * mass / 2.f);
}

float GibbsSampler::computeDeltaLL(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2)
{
    if (r1 == r2)
    {
        return gaps::algo::deltaLL(mDMatrix.nRow(), mDMatrix.colPtr(r1),
            mSMatrix.colPtr(r1), mAPMatrix.colPtr(r1), mOtherMatrix->colPtr(c1),
            m1, mOtherMatrix->colPtr(c2), m2);
    }
    return computeDeltaLL(r1, c1, m1) + computeDeltaLL(r2, c2, m2);
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

unsigned GibbsSampler::dataRows() const
{
    return mDMatrix.nRow();
}

unsigned GibbsSampler::dataCols() const
{
    return mDMatrix.nCol();
}

float GibbsSampler::chi2() const
{   
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}
  
uint64_t GibbsSampler::nAtoms() const
{   
    return mDomain.size();
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

float GibbsSampler::deathProb(uint64_t nAtoms) const
{
    double size = static_cast<double>(mDomainLength);
    double term1 = (size - static_cast<double>(nAtoms)) / size;
    double term2 = mAlpha * static_cast<double>(mNumBins) * term1;
    return static_cast<double>(nAtoms) / (static_cast<double>(nAtoms) + term2);
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

// add an atom at pos, calculate mass either with an exponential distribution
// or with the gibbs mass calculation
void GibbsSampler::birth()
{
    uint64_t pos = mDomain.randomFreePosition();
    AtomicProposal prop = AtomicProposal('B', pos);

    unsigned row = getRow(prop.birthPos);
    unsigned col = getCol(prop.birthPos);

    // calculate proposed mass
    float mass = 0.f;
    if (canUseGibbs(col))
    {
        AlphaParameters alpha = alphaParameters(row, col);
        mass = gibbsMass(alpha, &(prop.rng)).value; // 0 if it fails
    }
    else
    {
        mass = prop.rng.exponential(mLambda);
    }

    // accept mass as long as it's non-zero
    if (mass >= gaps::epsilon)
    {
        mDomain.insert(prop.birthPos, mass);
        mMatrix(row, col) += mass;
        updateAPMatrix(row, col, mass);
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass calculation
void GibbsSampler::death()
{
    Atom* a = mDomain.randomAtom();
    AtomicProposal prop = AtomicProposal('D', a);

    // calculate bin for this atom
    unsigned row = getRow(prop.atom1->pos);
    unsigned col = getCol(prop.atom1->pos);

    // kill off atom
    float newVal = gaps::max(mMatrix(row, col) - prop.atom1->mass, 0.f);
    updateAPMatrix(row, col, newVal - mMatrix(row, col));
    mMatrix(row, col) = newVal;

    // calculate rebirth mass
    float rebirthMass = prop.atom1->mass;
    AlphaParameters alpha = alphaParameters(row, col);
    if (canUseGibbs(col))
    {
        OptionalFloat gMass = gibbsMass(alpha, &(prop.rng));
        if (gMass.hasValue)
        {
            rebirthMass = gMass.value;
        }
    }

    // accept/reject rebirth
    float deltaLL = rebirthMass * (alpha.su - alpha.s * rebirthMass / 2.f);
    if (deltaLL * mAnnealingTemp >= std::log(prop.rng.uniform()))
    {
        prop.atom1->mass = rebirthMass;
        mMatrix(row, col) += rebirthMass;
        updateAPMatrix(row, col, rebirthMass);
    }
    else
    {
        mDomain.erase(prop.atom1->pos);
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
        AtomicProposal prop = AtomicProposal('M', hood.center, newLocation);

        float deltaLL = computeDeltaLL(r1, c1, -1.f * prop.atom1->mass,
            r2, c2, prop.atom1->mass);
        if (deltaLL * mAnnealingTemp > std::log(prop.rng.uniform()))
        {
            prop.atom1->pos = prop.moveDest;

            float newVal = gaps::max(mMatrix(r1, c1) - prop.atom1->mass, 0.f);
            updateAPMatrix(r1, c1, newVal - mMatrix(r1, c1));
            mMatrix(r1, c1) = newVal;

            mMatrix(r2, c2) += prop.atom1->mass;
            updateAPMatrix(r2, c2, prop.atom1->mass);
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
        AtomicProposal prop = AtomicProposal('E', a1, a2);

        float m1 = prop.atom1->mass;
        float m2 = prop.atom2->mass;

        if (canUseGibbs(c1, c2))
        {
            AlphaParameters alpha = alphaParameters(r1, c1, r2, c2);
            OptionalFloat gMass = gibbsMass(alpha, m1, m2, &(prop.rng));
            if (gMass.hasValue)
            {
                acceptExchange(&prop, gMass.value, r1, c1, r2, c2);
                return;
            }
        }

        float newMass = prop.rng.truncGammaUpper(m1 + m2, 2.f, 1.f / mLambda);

        float delta = m1 > m2 ? newMass - m1 : m2 - newMass; // change larger mass
        float pOldMass = 2.f * newMass > m1 + m2 ? gaps::max(m1, m2) : gaps::min(m1, m2);

        float pNew = gaps::d_gamma(newMass, 2.f, 1.f / mLambda);
        float pOld = gaps::d_gamma(pOldMass, 2.f, 1.f / mLambda);

        if (pOld == 0.f && pNew != 0.f) // special case
        {
            acceptExchange(&prop, delta, r1, c1, r2, c2);
            return;
        }

        float deltaLL = computeDeltaLL(r1, c1, delta, r2, c2, -delta);
        float priorLL = (pOld == 0.f) ? 1.f : pOld / pNew;
        float u = std::log(prop.rng.uniform() * priorLL);
        if (u < deltaLL * mAnnealingTemp)
        {
            acceptExchange(&prop, delta, r1, c1, r2, c2);
            return;
        }
    }
}

// helper function for acceptExchange, used to conditionally update the mass
// at a single position, deleting the atom if the new mass is too small
bool GibbsSampler::updateAtomMass(Atom *atom, float delta)
{
    if (atom->mass + delta < gaps::epsilon)
    {
        DEBUG_PING // want to know if this ever happens
        mDomain.erase(atom->pos);
        return false;
    }
    atom->mass += delta;
    return true;
}

// helper function for exchange step, updates the atomic domain, matrix, and
// A*P matrix
void GibbsSampler::acceptExchange(AtomicProposal *prop,
float d1, unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    float d2 = -1.f * d1;
    bool b1 = updateAtomMass(prop->atom1, d1);
    bool b2 = updateAtomMass(prop->atom2, d2);
    GAPS_ASSERT(b1 || b2);
    
    // delete entire atom if resize would make it too small
    if (!b1) { d1 = -1.f * prop->atom1->mass; }
    if (!b2) { d2 = -1.f * prop->atom2->mass; }

    // ensure matrix values don't go negative (truncation error at fault)
    float v1 = gaps::max(mMatrix(r1, c1) + d1, 0.f);
    updateAPMatrix(r1, c1, v1 - mMatrix(r1, c1));
    mMatrix(r1, c1) = v1;


    float v2 = gaps::max(mMatrix(r2, c2) + d2, 0.f);
    updateAPMatrix(r2, c2, v2 - mMatrix(r2, c2));
    mMatrix(r2, c2) = v2;
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

#ifdef GAPS_DEBUG
bool GibbsSampler::internallyConsistent()
{
    return true;
}
#endif // GAPS_DEBUG

Archive& operator<<(Archive &ar, GibbsSampler &s)
{
    ar << s.mMatrix << s.mAPMatrix << s.mDomain << s.mLambda
        << s.mMaxGibbsMass << s.mAnnealingTemp << s.mNumRows << s.mNumCols
        << s.mBinSize << s.mAvgQueue << s.mNumQueues;
    return ar;
}

Archive& operator>>(Archive &ar, GibbsSampler &s)
{
    ar >> s.mMatrix >> s.mAPMatrix >> s.mDomain >> s.mLambda
        >> s.mMaxGibbsMass >> s.mAnnealingTemp >> s.mNumRows >> s.mNumCols
        >> s.mBinSize >> s.mAvgQueue >> s.mNumQueues;
    return ar;
}