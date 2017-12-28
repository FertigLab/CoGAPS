#include "GibbsSampler.h"
#include "Algorithms.h"

#include <Rcpp.h>

#ifdef GAPS_DEBUG
#include <stdexcept>
#include <iostream>
#endif

static const double EPSILON = 1.e-10;
static const double AUTO_ACCEPT = std::numeric_limits<double>::min();

GibbsSampler::GibbsSampler(Rcpp::NumericMatrix D, Rcpp::NumericMatrix S,
unsigned int nFactor, double alphaA, double alphaP, double maxGibbsMassA,
double maxGibbsMassP, bool singleCellRNASeq)
    :
mDMatrix(D), mSMatrix(S), mAMatrix(D.nrow(), nFactor),
mPMatrix(nFactor, D.ncol()), mAPMatrix(D.nrow(), D.ncol()),
mADomain('A', D.nrow(), nFactor), mPDomain('P', nFactor, D.ncol()),
mMaxGibbsMassA(maxGibbsMassA), mMaxGibbsMassP(maxGibbsMassP),
mAnnealingTemp(1.0), mChi2(0.0), mSingleCellRNASeq(singleCellRNASeq),
mAMeanMatrix(D.nrow(), nFactor), mAStdMatrix(D.nrow(), nFactor),
mPMeanMatrix(nFactor, D.ncol()), mPStdMatrix(nFactor, D.ncol()),
mStatUpdates(0)
{
    double meanD = mSingleCellRNASeq ? gaps::algo::nonZeroMean(mDMatrix)
        : gaps::algo::mean(mDMatrix);

    mADomain.setAlpha(alphaA);
    mADomain.setLambda(alphaA * std::sqrt(nFactor / meanD));
    mPDomain.setAlpha(alphaP);
    mPDomain.setLambda(alphaP * std::sqrt(nFactor / meanD));

    mMaxGibbsMassA /= mADomain.lambda();
    mMaxGibbsMassP /= mPDomain.lambda();

    mChi2 = 2.0 * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}

double GibbsSampler::getGibbsMass(const MatrixChange &change)
{
    // check if this change is death (only called in birth/death)
    bool death = change.delta1 < 0;

    // get s and su
    AlphaParameters alphaParam = gaps::algo::alphaParameters(change, mDMatrix,
        mSMatrix, mAMatrix, mPMatrix, mAPMatrix);

    // calculate mean and standard deviation
    alphaParam.s *= mAnnealingTemp / 2.0;
    alphaParam.su *= mAnnealingTemp / 2.0;
    double lambda = change.label == 'A' ? mADomain.lambda() : mPDomain.lambda();
    double mean  = (2.0 * alphaParam.su - lambda) / (2.0 * alphaParam.s);
    double sd = 1.0 / std::sqrt(2.0 * alphaParam.s);

    // note: is bounded below by zero so have to use inverse sampling!
    // based upon algorithm in DistScalarRmath.cc (scalarRandomSample)
    double plower = gaps::random::p_norm(0.0, mean, sd);
    double u = gaps::random::uniform(plower, 1.0);

    // if the likelihood is flat and nonzero, sample strictly from the prior
    double newMass = 0.0;
    if (plower == 1 || alphaParam.s < 1.e-5)
    {
        newMass = death ? std::abs(change.delta1) : 0.0;
    }
    else if (plower >= 0.99)
    {
        double tmp1 = gaps::random::d_norm(0, mean, sd);
        double tmp2 = gaps::random::d_norm(10 * lambda, mean, sd);

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

    return std::max(newMass, 0.0);
}

double GibbsSampler::computeDeltaLL(const MatrixChange &change)
{
    return gaps::algo::deltaLL(change, mDMatrix, mSMatrix, mAMatrix,
        mPMatrix, mAPMatrix);
}

void GibbsSampler::update(char matrixLabel)
{
    std::cout << "temp: " << mAnnealingTemp << '\n';
    AtomicSupport &domain(matrixLabel == 'A' ? mADomain : mPDomain);

    // get proposal and convert to a matrix update
    AtomicProposal proposal = domain.makeProposal();

    std::cout << matrixLabel << " " << proposal.type << " " << proposal.pos1 << " " <<
        proposal.delta1 << " " << proposal.pos2 << " " << proposal.delta2 << '\n';

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

void GibbsSampler::setChi2(double chi2)
{
    mChi2 = chi2;
}

double GibbsSampler::chi2() const
{
    return mChi2;
}

void GibbsSampler::setAnnealingTemp(double temp)
{
    mAnnealingTemp = temp;
}

bool GibbsSampler::evaluateChange(AtomicSupport &domain,
const AtomicProposal &proposal, double threshold, bool accept)
{
    MatrixChange change = domain.getMatrixChange(proposal);
    double delLL = computeDeltaLL(change);
    if (accept || delLL * mAnnealingTemp >= threshold)
    {
        change = domain.acceptProposal(proposal);
        change.label == 'A' ? mAMatrix.update(change) : mPMatrix.update(change);
        updateAPMatrix(change);
        setChi2(chi2() - 2 * delLL);
        return true;
    }
    return false;
}

void GibbsSampler::updateAPMatrix_A(unsigned row, unsigned col, double delta)
{
    double newVal = 0.0;
    for (unsigned j = 0; j < mAPMatrix.nCol(); ++j)
    {
        newVal = mAPMatrix(row,j) + delta * mPMatrix(col,j);
        mAPMatrix.set(row, j, newVal);
    }
}

void GibbsSampler::updateAPMatrix_P(unsigned row, unsigned col, double delta)
{
    double newVal = 0.0;
    for (unsigned i = 0; i < mAPMatrix.nRow(); ++i)
    {
        newVal = mAPMatrix(i,col) + delta * mAMatrix(i,row);
        mAPMatrix.set(i, col, newVal);
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

        std::cout << "Gibbs Check: " << ch.label << " " << !(check1 && check2) << '\n';
        return !(check1 && check2);
    }
    std::cout << "Gibbs Check: " << ch.label << " " << !check1 << '\n';
    return !check1;
}

// accept automatically, try to rebirth
bool GibbsSampler::death(AtomicSupport &domain, AtomicProposal &proposal)
{
    // automaticallly accept death
    evaluateChange(domain, proposal, 0.0, true);

    // rebirth
    MatrixChange ch = domain.getMatrixChange(proposal);
    double newMass = canUseGibbs(ch) ? getGibbsMass(ch) : 0.0;
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

    uint64_t pos1 = proposal.delta1 > 0 ? proposal.pos1 : proposal.pos2;
    uint64_t pos2 = proposal.delta1 > 0 ? proposal.pos2 : proposal.pos1;

    double mass1 = domain.at(pos1);
    double mass2 = domain.at(pos2);

    double newMass1 = mass1 + std::max(proposal.delta1, proposal.delta2);
    double newMass2 = mass2 + std::min(proposal.delta1, proposal.delta2);

    unsigned row1 = change.delta1 > 0 ? change.row1 : change.row2;
    unsigned row2 = change.delta1 > 0 ? change.row2 : change.row1;

    unsigned col1 = change.delta1 > 0 ? change.col1 : change.col2;
    unsigned col2 = change.delta1 > 0 ? change.col2 : change.col1;

    change.row1 = row1;
    change.col1 = col1;
    change.row2 = row2;
    change.col2 = col2;
    change.delta1 = newMass1 - mass1;
    change.delta2 = newMass2 - mass2;

    if (canUseGibbs(change))
    {
        AlphaParameters alphaParam = gaps::algo::alphaParameters(change,
            mDMatrix, mSMatrix, mAMatrix, mPMatrix, mAPMatrix);
        alphaParam.s *= mAnnealingTemp;
        alphaParam.su *= mAnnealingTemp;

        if (alphaParam.s > EPSILON)
        {
            double mean = alphaParam.su / alphaParam.s;
            double sd = 1.0 / std::sqrt(alphaParam.s);

            double plower = gaps::random::p_norm(-mass1, mean, sd);
            double pupper = gaps::random::p_norm(mass2, mean, sd);
            double u = gaps::random::uniform(plower, pupper);

            if (!(plower >  0.95 || pupper < 0.05))
            {
                proposal.delta1 = gaps::random::q_norm(u, mean, sd);
                proposal.delta1 = std::max(proposal.delta1, -mass1);
                proposal.delta1 = std::min(proposal.delta1, mass2);
                proposal.pos1 = pos1;
                proposal.pos2 = pos2;
                proposal.delta2 = -proposal.delta1;

                return evaluateChange(domain, proposal, 0.0, true);
            }
        }
    }

    double pnewMass = mass1 > mass2 ? newMass1 : newMass2;
    double poldMass = newMass1 > newMass2 ? mass1 : mass2;

    double pnew = gaps::random::d_gamma(pnewMass, 2.0, 1.0 / domain.lambda());
    double pold = gaps::random::d_gamma(poldMass, 2.0, 1.0 / domain.lambda());

    if (pold != 0.0 || pnew == 0.0)
    {
        double priorLL = pnew == 0.0 ? 0.0 : log(pnew / pold);
        proposal.pos1 = pos1;
        proposal.delta1 = newMass1 - mass1;
        proposal.pos2 = pos2;
        proposal.delta2 = newMass2 - mass2;
        return evaluateChange(domain, proposal, std::log(gaps::random::uniform())
            - priorLL);
    }
    return false;
}

Rcpp::NumericMatrix GibbsSampler::AMeanRMatrix() const
{
    return gaps::algo::scalarDivision(mAMeanMatrix, mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GibbsSampler::AStdRMatrix() const
{
    return gaps::algo::computeStdDev(mAStdMatrix, mAMeanMatrix,
        mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GibbsSampler::PMeanRMatrix() const
{
    return gaps::algo::scalarDivision(mPMeanMatrix, mStatUpdates).rMatrix();
}

Rcpp::NumericMatrix GibbsSampler::PStdRMatrix() const
{
    return gaps::algo::computeStdDev(mPStdMatrix, mPMeanMatrix,
        mStatUpdates).rMatrix();
}

// TODO cache matrices used in these calculations
void GibbsSampler::updateStatistics()
{
    mStatUpdates++;

    Vector normVec(mPMatrix.nRow());
        
    for (unsigned r = 0; r < mPMatrix.nRow(); ++r)
    {
        normVec(r) = gaps::algo::sum(mPMatrix.getRow(r));
        if (normVec(r) == 0)
        {
            normVec(r) = 1.0;
        }
    }

    for (unsigned c = 0; c < mAMeanMatrix.nCol(); ++c)
    {
        mAMeanMatrix.getCol(c) += gaps::algo::scalarMultiple(mAMatrix.getCol(c),
            normVec(c));

        mAStdMatrix.getCol(c) += gaps::algo::squaredScalarMultiple(mAMatrix.getCol(c),
            normVec(c));
    }

    for (unsigned r = 0; r < mPMeanMatrix.nRow(); ++r)
    {
        mPMeanMatrix.getRow(r) += gaps::algo::scalarMultiple(mPMatrix.getRow(r),
            1.0 / normVec(r));

        mPStdMatrix.getRow(r) += gaps::algo::squaredScalarMultiple(mPMatrix.getRow(r),
            1.0 / normVec(r));
    }
}

#ifdef GAPS_DEBUG

void GibbsSampler::checkAtomicMatrixConsistency() const
{
    double AMass = gaps::algo::sum(mAMatrix);
    double PMass = gaps::algo::sum(mPMatrix);
    double Adiff = std::abs(AMass - mADomain.totalMass());
    double Pdiff = std::abs(PMass - mPDomain.totalMass());
    if (Adiff > 0.00001 || Pdiff > 0.00001)
    {
        std::cout << "A Mass - " << AMass << "," << mADomain.totalMass() << '\n';
        std::cout << "P Mass - " << PMass << "," << mPDomain.totalMass() << '\n';
        throw std::runtime_error("inconsistency with Atomic Domain");
    }
}

void GibbsSampler::checkAPMatrix() const
{
    if (!gaps::algo::checkAPMatrix(mAPMatrix, mAMatrix, mPMatrix))
    {
        throw std::runtime_error("APMatrix inconsistency");
    }
}

#endif