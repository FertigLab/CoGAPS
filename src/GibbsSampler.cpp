#include "GibbsSampler.h"
#include "Algorithms.h"

#include <Rcpp.h>

static const float EPSILON = 1.e-10;

GibbsSampler::GibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor)
    :
mDMatrix(D), mSMatrix(S), mAPMatrix(D.nrow(), D.ncol()),
mAMatrix(D.nrow(), nFactor), mPMatrix(nFactor, D.ncol()),
mADomain('A', D.nrow(), nFactor), mPDomain('P', nFactor, D.ncol()),
mAMeanMatrix(D.nrow(), nFactor), mAStdMatrix(D.nrow(), nFactor),
mPMeanMatrix(nFactor, D.ncol()), mPStdMatrix(nFactor, D.ncol()),
mPumpMatrix(D.nrow(), nFactor)
{}

GibbsSampler::GibbsSampler(const Rcpp::NumericMatrix &D,
const Rcpp::NumericMatrix &S, unsigned nFactor, float alphaA, float alphaP,
float maxGibbmassA, float maxGibbmassP, bool singleCellRNASeq,
char whichMatrixFixed, const Rcpp::NumericMatrix &FP, PumpThreshold pumpThreshold)
    :
mDMatrix(D), mSMatrix(S), mAPMatrix(D.nrow(), D.ncol()),
mAMatrix(D.nrow(), nFactor), mPMatrix(nFactor, D.ncol()),
mADomain('A', D.nrow(), nFactor), mPDomain('P', nFactor, D.ncol()),
mAMeanMatrix(D.nrow(), nFactor), mAStdMatrix(D.nrow(), nFactor),
mPMeanMatrix(nFactor, D.ncol()), mPStdMatrix(nFactor, D.ncol()),
mPumpMatrix(D.nrow(), nFactor), mPumpThreshold(pumpThreshold), mStatUpdates(0),
mPumpStatUpdates(0), mMaxGibbsMassA(maxGibbmassA), mMaxGibbsMassP(maxGibbmassP),
mAnnealingTemp(1.0), mSingleCellRNASeq(singleCellRNASeq), mNumFixedPatterns(0),
mFixedMat(whichMatrixFixed)
{
    float meanD = mSingleCellRNASeq ? gaps::algo::nonZeroMean(mDMatrix)
        : gaps::algo::mean(mDMatrix);

    mADomain.setAlpha(alphaA);
    mADomain.setLambda(alphaA * std::sqrt(nFactor / meanD));
    mPDomain.setAlpha(alphaP);
    mPDomain.setLambda(alphaP * std::sqrt(nFactor / meanD));

    mMaxGibbsMassA /= mADomain.lambda();
    mMaxGibbsMassP /= mPDomain.lambda();

    if (mFixedMat == 'A')
    {
        mNumFixedPatterns = FP.ncol();
        ColMatrix temp(FP);
        for (unsigned j = 0; j < mNumFixedPatterns; ++j)
        {
            mAMatrix.getCol(j) = temp.getCol(j) / gaps::algo::sum(temp.getCol(j));
        }
    }
    else if (mFixedMat == 'P')
    {
        mNumFixedPatterns = FP.nrow();
        RowMatrix temp(FP);
        for (unsigned i = 0; i < mNumFixedPatterns; ++i)
        {
            mPMatrix.getRow(i) = temp.getRow(i) / gaps::algo::sum(temp.getRow(i));
        }
    }
    gaps::algo::matrixMultiplication(mAPMatrix, mAMatrix, mPMatrix);
}

float GibbsSampler::getGibbsMass(const MatrixChange &change)
{
    // check if this change is death (only called in birth/death)
    bool death = change.delta1 < 0;

    // get s and su
    AlphaParameters alphaParam = gaps::algo::alphaParameters(change, mDMatrix,
        mSMatrix, mAMatrix, mPMatrix, mAPMatrix);

    // calculate mean and standard deviation
    alphaParam.s *= mAnnealingTemp / 2.0;
    alphaParam.su *= mAnnealingTemp / 2.0;
    float lambda = change.label == 'A' ? mADomain.lambda() : mPDomain.lambda();
    float mean  = (2.0 * alphaParam.su - lambda) / (2.0 * alphaParam.s);
    float sd = 1.0 / std::sqrt(2.0 * alphaParam.s);

    // note: is bounded below by zero so have to use inverse sampling!
    // based upon algorithm in DistScalarRmath.cc (scalarRandomSample)
    float plower = gaps::random::p_norm(0.f, mean, sd);
    float u = gaps::random::uniform(plower, 1.f);

    // if the likelihood is flat and nonzero, sample strictly from the prior
    float newMass = 0.f;
    if (plower == 1.f || alphaParam.s < 0.00001f)
    {
        newMass = death ? std::abs(change.delta1) : 0.f;
    }
    else if (plower >= 0.99f)
    {
        float tmp1 = gaps::random::d_norm(0.f, mean, sd);
        float tmp2 = gaps::random::d_norm(10.f * lambda, mean, sd);

        if (tmp1 > EPSILON && std::abs(tmp1 - tmp2) < EPSILON)
        {
            return death ? 0.0 : change.delta1;
        }
    }
    else
    {
        newMass = gaps::random::q_norm(u, mean, sd);
    }

    newMass = (change.label == 'A' ? std::min(newMass, mMaxGibbsMassA)
        : std::min(newMass, mMaxGibbsMassP));

    return std::max(newMass, 0.f);
}

float GibbsSampler::computeDeltaLL(const MatrixChange &change)
{
    return gaps::algo::deltaLL(change, mDMatrix, mSMatrix, mAMatrix,
        mPMatrix, mAPMatrix);
}

void GibbsSampler::update(char matrixLabel)
{
    AtomicSupport &domain(matrixLabel == 'A' ? mADomain : mPDomain);

    // get proposal and convert to a matrix update
    AtomicProposal proposal = domain.makeProposal();

    // check if proposal is in fixed region
    if (mFixedMat == 'A' || mFixedMat == 'P')
    {
        MatrixChange change = domain.getMatrixChange(proposal);
        unsigned index1 = mFixedMat == 'A' ? change.col1 : change.row1;
        unsigned index2 = mFixedMat == 'A' ? change.col2 : change.row2;
        index2 = change.nChanges == 1 ? mNumFixedPatterns : index2;

        if (index1 < mNumFixedPatterns || index2 < mNumFixedPatterns)
        {
            return;
        }
    }

    // Update based on the proposal type
    switch (proposal.type)
    {
        case 'D': death(domain, proposal);    break;
        case 'B': birth(domain, proposal);    break;
        case 'M': move(domain, proposal);     break;
        case 'E': exchange(domain, proposal); break;
    }
}

uint64_t GibbsSampler::totalNumAtoms(char matrixLabel) const
{
    return matrixLabel == 'A' ? mADomain.numAtoms() : mPDomain.numAtoms();
}

float GibbsSampler::chi2() const
{
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}

void GibbsSampler::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

bool GibbsSampler::evaluateChange(AtomicSupport &domain,
const AtomicProposal &proposal, float threshold, bool accept)
{
    MatrixChange change = domain.getMatrixChange(proposal);
    float delLL = accept ? 0.f : computeDeltaLL(change);
    if (accept || delLL * mAnnealingTemp >= threshold)
    {
        change = domain.acceptProposal(proposal);
        change.label == 'A' ? mAMatrix.update(change) : mPMatrix.update(change);
        updateAPMatrix(change);
        return true;
    }
    return false;
}

// simd?
void GibbsSampler::updateAPMatrix_A(unsigned row, unsigned col, float delta)
{
    const Vector &APvec = mAPMatrix.getRow(row);
    const Vector &Pvec = mPMatrix.getRow(col);
    for (unsigned j = 0; j < mAPMatrix.nCol(); ++j)
    {
        mAPMatrix.set(row, j, APvec[j] + delta * Pvec[j]);
    }
}

// simd?
void GibbsSampler::updateAPMatrix_P(unsigned row, unsigned col, float delta)
{
    const Vector &APvec = mAPMatrix.getCol(col);
    const Vector &Avec = mAMatrix.getCol(row);
    for (unsigned i = 0; i < mAPMatrix.nRow(); ++i)
    {
        mAPMatrix.set(i, col, APvec[i] + delta * Avec[i]);
    }
}

void GibbsSampler::updateAPMatrix(const MatrixChange &change)
{
    if (change.label == 'A')
    {
        updateAPMatrix_A(change.row1, change.col1, change.delta1);
        if (change.nChanges > 1)
        {
            updateAPMatrix_A(change.row2, change.col2, change.delta2);
        }
    }
    else
    {
        updateAPMatrix_P(change.row1, change.col1, change.delta1);
        if (change.nChanges > 1)
        {
            updateAPMatrix_P(change.row2, change.col2, change.delta2);
        }
    }
}

bool GibbsSampler::canUseGibbs(const MatrixChange &ch)
{
    bool check1 = (ch.label == 'A' && gaps::algo::isRowZero(mPMatrix, ch.col1))
        || (ch.label == 'P' && gaps::algo::isColZero(mAMatrix, ch.row1));

    if (ch.nChanges > 1)
    {
        bool check2 = (ch.label == 'A' && gaps::algo::isRowZero(mPMatrix, ch.col2))
            || (ch.label == 'P' && gaps::algo::isColZero(mAMatrix, ch.row2));
        return !(check1 && check2);
    }
    return !check1;
}

// accept automatically, try to rebirth
bool GibbsSampler::death(AtomicSupport &domain, AtomicProposal &proposal)
{
    // automaticallly accept death
    evaluateChange(domain, proposal, 0.0, true);

    // rebirth
    MatrixChange ch = domain.getMatrixChange(proposal);
    float newMass = canUseGibbs(ch) ? getGibbsMass(ch) : 0.0;
    proposal.delta1 = newMass < EPSILON ? -proposal.delta1 : newMass;    

    // attempt to accept rebirth
    return evaluateChange(domain, proposal, std::log(gaps::random::uniform()));
}

bool GibbsSampler::birth(AtomicSupport &domain, AtomicProposal &proposal)
{
    // attempt gibbs
    MatrixChange ch = domain.getMatrixChange(proposal);
    proposal.delta1 = canUseGibbs(ch) ? getGibbsMass(ch) : proposal.delta1;

    // accept birth
    return evaluateChange(domain, proposal, 0.0, true);
}

bool GibbsSampler::move(AtomicSupport &domain, AtomicProposal &proposal)
{
    MatrixChange change = domain.getMatrixChange(proposal);
    if (change.row1 == change.row2 && change.col1 == change.col2)
    {
        return false;
    }
    return evaluateChange(domain, proposal, std::log(gaps::random::uniform()));
}

bool GibbsSampler::exchange(AtomicSupport &domain, AtomicProposal &proposal)
{
    MatrixChange change = domain.getMatrixChange(proposal);
    if (change.row1 == change.row2 && change.col1 == change.col2)
    {
        return false;
    }

    float mass1 = domain.at(proposal.pos1);
    float mass2 = domain.at(proposal.pos2);
    float newMass1 = mass1 + proposal.delta1;
    float newMass2 = mass2 + proposal.delta2;

    if (canUseGibbs(change))
    {
        AlphaParameters alphaParam = gaps::algo::alphaParameters(change,
            mDMatrix, mSMatrix, mAMatrix, mPMatrix, mAPMatrix);
        alphaParam.s *= mAnnealingTemp;
        alphaParam.su *= mAnnealingTemp;

        if (alphaParam.s > EPSILON)
        {
            float mean = alphaParam.su / alphaParam.s;
            float sd = 1.f / std::sqrt(alphaParam.s);

            float plower = gaps::random::p_norm(-mass1, mean, sd);
            float pupper = gaps::random::p_norm(mass2, mean, sd);
            float u = gaps::random::uniform(plower, pupper);

            if (!(plower >  0.95f || pupper < 0.05f))
            {
                proposal.delta1 = gaps::random::q_norm(u, mean, sd);
                proposal.delta1 = std::max(proposal.delta1, -mass1);
                proposal.delta1 = std::min(proposal.delta1, mass2);
                proposal.delta2 = -proposal.delta1;
                return evaluateChange(domain, proposal, 0.0, true);
            }
        }
    }

    float pnewMass = mass1 > mass2 ? newMass1 : newMass2;
    float poldMass = newMass1 > newMass2 ? mass1 : mass2;

    float pnew = gaps::random::d_gamma(pnewMass, 2.0, 1.f / domain.lambda());
    float pold = gaps::random::d_gamma(poldMass, 2.0, 1.f / domain.lambda());

    if (pold == 0.f && pnew != 0.f)
    {
        return evaluateChange(domain, proposal, 0.f, true);
    }
    else
    {
        float priorLL = (pold == 0.f) ? 0.f : log(pnew / pold);
        float rejectProb = std::log(gaps::random::uniform()) - priorLL;
        return evaluateChange(domain, proposal, rejectProb);
    }
}

Rcpp::NumericMatrix GibbsSampler::AMeanRMatrix() const
{
    return (mAMeanMatrix / mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GibbsSampler::AStdRMatrix() const
{
    return gaps::algo::computeStdDev(mAStdMatrix, mAMeanMatrix,
        mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GibbsSampler::PMeanRMatrix() const
{
    return (mPMeanMatrix / mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GibbsSampler::PStdRMatrix() const
{
    return gaps::algo::computeStdDev(mPStdMatrix, mPMeanMatrix,
        mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GibbsSampler::pumpMatrix() const
{
    unsigned denom = mPumpStatUpdates ? mPumpStatUpdates : 1.f;
    return (mPumpMatrix / denom).rMatrix();
}

Rcpp::NumericMatrix GibbsSampler::meanPattern()
{
    ColMatrix Amean(mAMeanMatrix / (float)mStatUpdates);
    RowMatrix Pmean(mPMeanMatrix / (float)mStatUpdates);
    ColMatrix mat(mAMatrix.nRow(), mAMatrix.nCol());
    patternMarkers(Amean, Pmean, mat);
    return mat.rMatrix();
}

float GibbsSampler::meanChiSq() const
{
    ColMatrix Amean = mAMeanMatrix / (float)mStatUpdates;
    RowMatrix Pmean = mPMeanMatrix / (float)mStatUpdates;
    TwoWayMatrix Mmean(Amean.nRow(), Pmean.nCol());
    gaps::algo::matrixMultiplication(Mmean, Amean, Pmean);
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, Mmean);
}

void GibbsSampler::updateStatistics()
{
    mStatUpdates++;
    unsigned nPatterns = mAMatrix.nCol();

    Vector normVec(nPatterns);
    for (unsigned j = 0; j < nPatterns; ++j)
    {
        normVec[j] = gaps::algo::sum(mPMatrix.getRow(j));
        normVec[j] = normVec[j] == 0 ? 1.f : normVec[j];

        Vector quot(mPMatrix.getRow(j) / normVec[j]);
        mPMeanMatrix.getRow(j) += quot;
        mPStdMatrix.getRow(j) += gaps::algo::elementSq(quot);

        Vector prod(mAMatrix.getCol(j) * normVec[j]);
        mAMeanMatrix.getCol(j) += prod;
        mAStdMatrix.getCol(j) += gaps::algo::elementSq(prod); 
    }
}

void GibbsSampler::updatePumpStatistics()
{
    if (mFixedMat != 'A')
    {
        mPumpStatUpdates++;
        patternMarkers(normedAMatrix(), normedPMatrix(), mPumpMatrix);
    }
}

ColMatrix GibbsSampler::normedAMatrix() const
{
    ColMatrix normedA(mAMatrix);
    for (unsigned j = 0; j < normedA.nCol(); ++j)
    {
        float factor = gaps::algo::sum(mPMatrix.getRow(j));
        factor = (factor == 0) ? 1.f : factor;
        normedA.getCol(j) = normedA.getCol(j) * factor;
    }
    return normedA;
}

RowMatrix GibbsSampler::normedPMatrix() const
{
    RowMatrix normedP(mPMatrix);
    for (unsigned i = 0; i < normedP.nRow(); ++i)
    {
        float factor = gaps::algo::sum(mPMatrix.getRow(i));
        factor = (factor == 0) ? 1.f : factor;
        normedP.getRow(i) = normedP.getRow(i) / factor;
    }
    return normedP;
}

static unsigned geneThreshold(const ColMatrix &rankMatrix, unsigned pat)
{
    float cutRank = rankMatrix.nRow();
    for (unsigned i = 0; i < rankMatrix.nRow(); ++i)
    {
        for (unsigned j = 0; j < rankMatrix.nCol(); ++j)
        {
            if (j != pat && rankMatrix(i,j) <= rankMatrix(i,pat))
            {
                cutRank = std::min(cutRank, std::max(0.f, rankMatrix(i,pat)-1));
            }
        }
    }
    return static_cast<unsigned>(cutRank);
}

void GibbsSampler::patternMarkers(RowMatrix normedA, RowMatrix normedP,
ColMatrix &statMatrix)
{
    // helpful notation
    unsigned nGenes = normedA.nRow();
    unsigned nPatterns = normedA.nCol();

    // scale A matrix
    for (unsigned j = 0; j < nPatterns; ++j)
    {
        float scale = gaps::algo::max(normedP.getRow(j));
        for (unsigned i = 0; i < nGenes; ++i)
        {
            normedA(i,j) *= scale;
        }
    }
    
    // compute sstat
    TwoWayMatrix sStat(nGenes, nPatterns);
    Vector lp(nPatterns), diff(nPatterns);
    for (unsigned j = 0; j < nPatterns; ++j)
    {
        lp[j] = 1.f;
        for (unsigned i = 0; i < nGenes; ++i)
        {
            float geneMax = gaps::algo::max(normedA.getRow(i));
            diff = geneMax > 0.f ? normedA.getRow(i) / geneMax - lp : lp * -1.f;
            sStat.set(i, j, std::sqrt(gaps::algo::dot(diff, diff)));
        }
        lp[j] = 0.f;
    }

    // update PUMP matrix
    if (mPumpThreshold == PUMP_UNIQUE)
    {
        for (unsigned i = 0; i < nGenes; ++i)
        {
            unsigned minNdx = gaps::algo::whichMin(sStat.getRow(i));
            statMatrix(i,minNdx)++;
        }
    }
    else if (mPumpThreshold == PUMP_CUT)
    {
        ColMatrix rankMatrix(nGenes, nPatterns);
        for (unsigned j = 0; j < nPatterns; ++j)
        {
            rankMatrix.getCol(j) = gaps::algo::rank(sStat.getCol(j));
        }
        
        for (unsigned j = 0; j < nPatterns; ++j)
        {
            unsigned cutRank = geneThreshold(rankMatrix, j);
            for (unsigned i = 0; i < nGenes; ++i)
            {
                if (rankMatrix(i,j) <= cutRank)
                {
                    statMatrix(i,j)++;
                }
            }
        }
    }
}

Archive& operator<<(Archive &ar, GibbsSampler &sampler)
{
    ar << sampler.mAMatrix << sampler.mPMatrix << sampler.mADomain
        << sampler.mPDomain << sampler.mAMeanMatrix << sampler.mAStdMatrix
        << sampler.mPMeanMatrix << sampler.mPStdMatrix << sampler.mStatUpdates
        << sampler.mPumpMatrix << sampler.mPumpThreshold << sampler.mPumpStatUpdates
        << sampler.mMaxGibbsMassA << sampler.mMaxGibbsMassP
        << sampler.mAnnealingTemp << sampler.mSingleCellRNASeq
        << sampler.mNumFixedPatterns << sampler.mFixedMat;

    return ar;
}

Archive& operator>>(Archive &ar, GibbsSampler &sampler)
{
    ar >> sampler.mAMatrix >> sampler.mPMatrix >> sampler.mADomain
        >> sampler.mPDomain >> sampler.mAMeanMatrix >> sampler.mAStdMatrix
        >> sampler.mPMeanMatrix >> sampler.mPStdMatrix >> sampler.mStatUpdates
        >> sampler.mPumpMatrix >> sampler.mPumpThreshold >> sampler.mPumpStatUpdates
        >> sampler.mMaxGibbsMassA >> sampler.mMaxGibbsMassP
        >> sampler.mAnnealingTemp >> sampler.mSingleCellRNASeq
        >> sampler.mNumFixedPatterns >> sampler.mFixedMat;

    gaps::algo::matrixMultiplication(sampler.mAPMatrix, sampler.mAMatrix,
        sampler.mPMatrix);
    
    return ar;
}