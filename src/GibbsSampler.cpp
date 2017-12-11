#include "GibbsSampler.h"
#include "Algorithms.h"

#include <Rcpp.h>

#define EPSILON 1E-10

GibbsSampler::GibbsSampler(Rcpp::NumericMatrix D, Rcpp::NumericMatrix S,
unsigned int nFactor, double alphaA, double alphaP, double maxGibbsMassA,
double maxGibbsMassP, bool singleCellRNASeq)
    :
mDMatrix(D), mSMatrix(S), mAMatrix(D.nrow(), nFactor),
mPMatrix(nFactor, D.ncol()), mAPMatrix(D.nrow(), D.ncol()),
mADomain('A', D.nrow(), nFactor), mPDomain('P', nFactor, D.ncol()),
mMaxGibbsMassA(maxGibbsMassA), mMaxGibbsMassP(maxGibbsMassP),
mAnnealingTemp(0.0), mChi2(0.0), mSingleCellRNASeq(singleCellRNASeq)
{
    double meanD = mSingleCellRNASeq ? gaps::algo::nonZeroMean(mDMatrix)
        : gaps::algo::mean(mDMatrix);

    mADomain.setAlpha(alphaA);
    mADomain.setLambda(alphaA * sqrt(nFactor / meanD));
    mPDomain.setAlpha(alphaP);
    mPDomain.setLambda(alphaP * sqrt(nFactor / meanD));
}

double GibbsSampler::getGibbsMass(const MatrixChange &change)
{
    // get s and su
    AlphaParameters alphaParam = gaps::algo::alphaParameters(change, mDMatrix,
        mSMatrix, mAMatrix, mPMatrix, mAPMatrix);

    // calculate mean and standard deviation
    double lambda = change.label == 'A' ? mADomain.lambda() : mPDomain.lambda();
    double mean  = (2.0 * alphaParam.su - lambda) / (2.0 * alphaParam.s);
    double sd = 1.0 / sqrt(2 * alphaParam.s);

    // note: is bounded below by zero so have to use inverse sampling!
    // based upon algorithm in DistScalarRmath.cc (scalarRandomSample)
    double plower = gaps::random::p_norm(0.0, mean, sd);
    double u = gaps::random::uniform(plower, 1.0);

    // if the likelihood is flat and nonzero, sample strictly from the prior
    double newMass = 0.0;
    if (plower == 1 || alphaParam.s < 1.e-5)
    {
        newMass = change.delta1 < 0 ? abs(change.delta1) : 0.0;
    }
    else if (plower >= 0.99)
    {
        double tmp1 = gaps::random::d_norm(0, mean, sd);
        double tmp2 = gaps::random::d_norm(10 * lambda, mean, sd);

        if (tmp1 > EPSILON && fabs(tmp1 - tmp2) < EPSILON)
        {
            return change.delta1 < 0 ? 0.0 : change.delta1;
        }
    }
    else
    {
        newMass = gaps::random::q_norm(u, mean, sd);
    }

    newMass = change.label == 'A' ? std::min(newMass, mMaxGibbsMassA)
        : std::min(newMass, mMaxGibbsMassP);

    // TODO check if less than epsilon and return original mass
    return std::max(newMass, 0.0);
}

double GibbsSampler::computeDeltaLL(const MatrixChange &change)
{
    return gaps::algo::deltaLL(change, mDMatrix, mSMatrix, mAMatrix,
        mPMatrix, mAPMatrix);
}

void GibbsSampler::update(char matrixLabel)
{
    AtomicSupport &domain(matrixLabel == 'A' ? mADomain : mPDomain);

    // get proposal and convert to a matrix update
    AtomicProposal proposal = domain.makeProposal();

    // Update based on the proposal type
    switch (proposal.type)
    {
        case 'D':
            death(domain, proposal);
            break;
        case 'B':
            birth(domain, proposal);
            break;
        case 'M':
            move(domain, proposal);
            break;
        case 'E':
            exchange(domain, proposal);
            break;
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
const AtomicProposal &proposal, double prob)
{
    MatrixChange change = domain.getMatrixChange(proposal);
    double delLL = computeDeltaLL(change);
    if (1 - delLL * mAnnealingTemp < prob)
    {
        domain.acceptProposal(proposal);
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

        return !check1 && !check2;
    }
    return !check1;
}

// accept automatically, try to rebirth
bool GibbsSampler::death(AtomicSupport &domain, AtomicProposal &proposal)
{
    // automaticallly accept death
    evaluateChange(domain, proposal, 1.0);

    // rebirth
    MatrixChange ch = domain.getMatrixChange(proposal);
    proposal.delta1 = canUseGibbs(ch) ? getGibbsMass(ch) : -proposal.delta1;

    // attempt to accept rebirth
    return evaluateChange(domain, proposal, log(gaps::random::uniform()));
}

bool GibbsSampler::birth(AtomicSupport &domain, AtomicProposal &proposal)
{
    // attempt gibbs
    MatrixChange change = domain.getMatrixChange(proposal);
    if (canUseGibbs(change))
    {
        proposal.delta1 = getGibbsMass(change);
    }
    return evaluateChange(domain, proposal, 1.0);
}

bool GibbsSampler::move(AtomicSupport &domain, AtomicProposal &proposal)
{
    MatrixChange change = domain.getMatrixChange(proposal);
    if (change.row1 == change.row2 && change.col1 == change.col2)
    {
        return false;
    }
    return evaluateChange(domain, proposal, log(gaps::random::uniform()));
}

bool GibbsSampler::exchange(AtomicSupport &domain, AtomicProposal &proposal)
{
    MatrixChange change = domain.getMatrixChange(proposal);
    if (change.row1 == change.row2 && change.col1 == change.col2)
    {
        return false;
    }


/*




    bool useGibbs = the_matrix_label == 'A' ? checkOtherMatrix('A', iRow, iCol, POrig) : checkOtherMatrix('P', iRow, iCol, AOrig);

    // -------------------------------------------------------------------------
    // EXCHANGE ACTION when initial useGibbs = true and check if Gibbs is usable.
    double s = 0, su = 0, mean = 0, sd = 0; // EJF -- MFO check
    std::pair<double, double> alphaparam(0,0);

    if (useGibbs)
    {
        alphaparam = GAPSNorm::calcAlphaParameters(the_matrix_label, _nFactor, D, S, AOrig, POrig, iGene1,
            iPattern1, iGene2, iPattern2, iSample1, iSample2);
    }

    if (alphaparam.first == 0.0)
    {
        useGibbs = false;
    }

    if (useGibbs)
    {
        s = alphaparam.first * _annealingTemperature;
        su = alphaparam.second * _annealingTemperature;
        mean = su / s;
        sd = 1. / sqrt(s);

        // set newMass1
        // need to retain exponential prior
        double plower = gaps::random::p_norm(-mass1, mean, sd);
        double pupper = gaps::random::p_norm(mass2, mean, sd);
        double u = plower + gaps::random::uniform() * (pupper - plower);

        // must sample from prior if the computed parameters are not good for Gibbs
        if (!(plower >  0.95 || pupper < 0.05 || s < epsilon ||
        newMass1 == DOUBLE_POSINF || newMass1 == DOUBLE_NEGINF))
        {
            proposal.deltaMass1 = gaps::random::q_norm(u, mean, sd);
            proposal.deltaMass1 = std::max(proposal.deltaMass1, -mass1);
            proposal.deltaMass1 = std::min(proposal.deltaMass1, mass2);
            proposal.deltaMass2 = -proposal.deltaMass2;

            std::vector<ElementChange> matChange = domain.getElementChange(proposal);
            double delLLnew = computeDeltaLL(the_matrix_label, D, S, AOrig, POrig, matChange);
            domain.acceptProposal(proposal);
            update_sysChi2(delLLnew);
            return true;
        }
    }

    // We can't use Gibbs, need
    // Metropolis-Hasting exchange action
    
    double lambda = the_matrix_label == 'A' ? _lambdaA : _lambdaP;
    double pnew = mass1 > mass2 ? gaps::random::d_gamma(newMass1, 2.0, 1.0 / lambda)
        : gaps::random::d_gamma(newMass2, 2.0, 1.0 / lambda)
    
    double pold = newMass1 > newMass2 ? gaps::random::d_gamma(mass1, 2.0, 1.0 / lambda)
        : gaps::random::d_gamma(mass2, 2.0, 1.0 / lambda);

    if (pnew == 0.0 && pold == 0.0)
    {
        priorLL = 0.0;
    }
    else if (pnew != 0.0 && pold == 0.0)
    {
        priorLL = DOUBLE_POSINF;
    }
    else
    {
        priorLL = log(pnew / pold);
    }

    double delLLnew = computeDeltaLL(the_matrix_label, D, S, AOrig, POrig, _matrixElemChange);

    proposal.deltaMass1 = newMass1 - mass1;
    proposal.deltaMass2 = newMass2 - mass2;
    std::vector<ElementChange> matChange = domain.getElementChange(proposal);
    
    if (priorLL == DOUBLE_POSINF) // update chi2 ??
    {
        return true;
    }

    if (1 - priorLL + delLLnew * _annealingTemperature  < log(gaps::random::uniform()))
    {
        domain.acceptProposal(false);
        update_sysChi2(delLLnew);
        return true;
    }
    return false;*/
}