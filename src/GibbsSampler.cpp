#include "GibbsSampler.h"

#include <Rcpp.h>

#if 0

GibbsSampler::GibbsSampler(Rcpp::NumericMatrix D, Rcpp::NumericMatrix S,
unsigned int nFactor, double alphaA, double alphaP)
    :
mDMatrix(D), mSMatrix(S), mAMatrix(D.nrow(), nFactor),
mPMatrix(nFactor, D.ncol()), mAPMatrix(D.nrow(), D.ncol()),
mADomain(D.nrow(), nFactor), mPDomain(nFactor, D.ncol())
{
    double meanD = mSingleCellRNASeq ? gaps::mat_algo::nonZeroMean(mDMatrix)
        : gaps::mat_algo::mean(mDMatrix);

    mADomain.setAlpha(alphaA);
    mADomain.setLambda(alphaA * sqrt(nFactor / meanD) * lambdaScaleFactorA);
    mPDomain.setAlpha(alphaP);
    mPDomain.setLambda(alphaP * sqrt(nFactor / meanD) * lambdaScaleFactorP);
}

double GibbsSampler::getGibbsMass(char matrixLabel, double origMass,
unsigned int iRow, unsigned int iCol, const Matrix &otherMatrix,
const Matrix &currentChainMatrix, const Matrix &D, const Matrix &S,
double rng)
{
    double s  = 0.0;
    double su = 0.0;

    if (matrixLabel == 'A')
    {
        double psq = 0.0, ssq = 0.0;
        for (unsigned j = 0; j < mDMatrix.nCol(); ++j)
        {
            psq = mPMatrix(col,j) * mPMatrix(col,j);
            ssq = mSMatrix(row,j) * mSMatrix(row,j);
           
            s += (mAnnealingTemp * psq / ssq) / 2;
            su += s * (mDMatrix(row,j) - mAPMatrix(row,j));
        }
    }
    else
    {
        double asq = 0.0, ssq = 0.0;
        for (unsigned i = 0; i < mDMatrix.nRow(); ++i)
        {
            asq = mAMatrix(i,row) * mAMatrix(i,row);
            ssq = mSMatrix(i,col) * mSMatrix(i,col);
           
            s += (mAnnealingTemp * asq / ssq) / 2;
            su += s * (mDMatrix(i,col) - mAPMatrix(i,col));
        }
    }

    // calculate mean and standard deviation
    double mean  = (2.0 * su - lambda) / (2.0 * s);
    double sd = 1.0 / sqrt(2 * s);

    // note: is bounded below by zero so have to use inverse sampling!
    // based upon algorithm in DistScalarRmath.cc (scalarRandomSample)
    double plower = gaps::random::p_norm(0.0, mean, sd);
    double u = gaps::random::uniform(plower, 1.0);

    // if the likelihood is flat and nonzero,
    // force to sample strictly from the prior
    if (plower == 1 || s < 1.e-5)
    {
        newMass = origMass < 0 ? abs(origMass) : 0.0;
    }
    else if (plower >= 0.99)
    {
        double tmp1 = Random::d_norm(0, mean, sd);
        double tmp2 = Random::d_norm(10 * lambda, mean, sd);

        if ((tmp1 > epsilon) && (fabs(tmp1 - tmp2) < epsilon))
        {
            return origMass < 0 ? 0.0 : origMass;
        }
    }
    else
    {
        newMass = Random::q_norm(u, mean, sd);
    }

    if (matrixLabel == 'A')
        newMass = std::min(newMass, mMaxGibbsMassA);
    else
        newMass = std::min(newMass mMaxGibbsMassP);

    return std::max(newMass, 0);
}

double GibbsSampler::cal_logLikelihood()
{
    return GAPSNorm::calChi2(mDMatrix, mAPMatrix, mSMatrix) / 2.0;
}

// -----------------------------------------------------------------------------
double GibbsSampler::computeDeltaLL(char the_matrix_label, const Matrix &D,
const Matrix &S, const Matrix &A, const Matrix &P,
const std::vector<ElementChange> the_matrixElemChange)
{
    if (the_matrixElemChange.size() == 0)
    {
        return 0.0;
    }
    else if (the_matrixElemChange.size() == 1)
    {
        return GAPSNorm::calcDeltaLL1E(the_matrix_label, D, S, A, P,
            the_matrixElemChange, _nFactor);
    }
    else if (the_matrixElemChange.size() == 2)
    {
        return GAPSNorm::calcDeltaLL2E(the_matrix_label, D, S, A, P,
            the_matrixElemChange, _nFactor);
    }
    else
    {
        return GAPSNorm::calcDeltaLLGen(the_matrix_label, D, S, A, P,
            the_matrixElemChange, _nFactor);
    }
}

// -----------------------------------------------------------------------------
void GibbsSampler::update(char label)
{
    const AtomicSupport &domain(label == 'A' ? _AAtomicdomain : _PAtomicdomain);
    const Matrix &matrix(label == 'A' ? _AMatrix : _PMatrix)

    // get proposal and convert to a matrix update
    AtomicProposal proposal = domain.makeProposal();
    std::vector<ElementChange> matChange = domain.getElementChange(proposal);

    // Update based on the proposal type
    bool Q_update;
    switch (proposal.type)
    {
        case 'D':
            Q_update = death(label, _DMatrix, _SMatrix, _AMatrix, _PMatrix, proposal, matChange);
            break;
        case 'B':
            Q_update = birth(label, _DMatrix, _SMatrix, _AMatrix, _PMatrix, proposal, matChange);
            break;
        case 'M':
            Q_update = move(label, _DMatrix, _SMatrix, _AMatrix, _PMatrix);
            break;
        case 'E':
            Q_update = exchange(label, _DMatrix, _SMatrix, _AMatrix, _PMatrix);
            break;
    }
    
    if (Q_update)
    {
        matrix.elemUpdate(matChange);
    }
}

// ----------------------------------------------------------------------------
void GibbsSampler::init_sysChi2()
{
    _sysChi2 = 2.0 * cal_logLikelihood();
}

void GibbsSampler::update_sysChi2(double delsysChi2)
{
    _sysChi2 -= 2.0 * delsysChi2;
}

double GibbsSampler::get_sysChi2()
{
    return _sysChi2;
}

// accept automatically, try to rebirth (effectively reject by attempting birth)
bool GibbsSampler::death(char the_matrix_label, const Matrix &D, const Matrix &S,
const Matrix &AOrig, const Matrix &POrig, const AtomicProposal &proposal,
const AtomicSupport &domain, const Matrix &matrix)
{
    std::vector<ElementChange> matChange = domain.getElementChange(proposal);

    unsigned row = matrixChange[0].row;
    unsigned col = matrixChange[0].col;

    double delLL = computeDeltaLL(matrixLabel, matrixChange);
    domain.acceptProposal(proposal);
    matrix.elemUpdate(matrixChange);
    update_sysChi2(delLL);

    // rebirth
    double newMass = -proposal.deltaMass1;
    bool gibbs = matrixLabel == 'A' ? gaps::mat_algo::isRowNonZero(mPMatrix, col)
        : gaps::mat_algo::isColNonZero(mAMatrix, row);

    if (gibbs)
    {
        newMass = the_matrix_label == 'A' ? getMass('A', origMass, iRow, iCol, POrig, AOrig, D, S, 0.1)
            : getMass('P', origMass, iRow, iCol, AOrig, POrig, D, S, 0.1);
        if (newMass <= epsilon)
        {
            newMass = -proposal.deltaMass1;
        }
    }

    proposal.deltaMass1 = newMass;
    std::vector<ElementChange> matChange = domain.getElementChange(proposal);                        
    delLL = computeDeltaLL('A', D, S, AOrig, POrig, matChange);

    // M-H sampling to determine whether or not we can accept Gibbs
    if (1 - delLLnew * _annealingTemperature  < log(Random::uniform()))
    {
        domain.acceptProposal(proposal);
        update_sysChi2(delLLnew);
        return true;
    }
    return false;
}

bool GibbsSampler::birth(char matrixLabel, const AtomicProposal &proposal,
const AtomicSupport &domain)
{
    double origMass = proposal.deltaMass1;
    std::vector<ElementChange> matChange = domain.getElementChange(proposal);

    unsigned int row = matrix[0].row;
    unsigned int col = matrix[0].col;

    bool gibbs = matrixLabel == 'A' ? canUseGibbs('A', iRow, iCol)
        : canUseGibbs('P', iRow, iCol);

    if (gibbs)
    {
        proposal.deltaMass1 = the_matrix_label == 'A' ? getMass('A', origMass, iRow, iCol, POrig, AOrig, D, S, 0.1)
            : getMass('P', origMass, iRow, iCol, AOrig, POrig, D, S, 0.1);
        std::vector<ElementChange> matChange = domain.getElementChange(proposal);                        
    }
    domain.acceptProposal(proposal);
    update_sysChi2(computeDeltaLL(the_matrix_label, D, S, AOrig, POrig, matChange)); 
    return true;
}

bool GibbsSampler::move(char the_matrix_label, const Matrix &D, const Matrix &S,
const Matrix &AOrig, const Matrix &POrig)
{
    std::vector<ElementChange> matChange = domain.getElementChange(proposal);

    if (matChange.row[0] == matChange[1].row && matChange[1].col == matChange[1].col)
    {
        return false;
    }

    double delLL = computeDeltaLL(matChange, D, S, AOrig, POrig, matChange);
    
    if (1 - delLL * _annealingTemperature < log(Random::uniform()))
    {
        domain.acceptProposal(proposal);
        update_sysChi2(delLL);
        return true;
    }
    return false;
}

bool GibbsSampler::exchange(char the_matrix_label, const Matrix &D,
const Matrix &S, const Matrix &AOrig, const Matrix &POrig)
{
    std::vector<ElementChange> matChange = domain.getElementChange(proposal);

    // return false if bin1 == bin2
    if (matChange.row[0] == matChange[1].row && matChange[1].col == matChange[1].col)
    {
        return false;
    }

    bool useGibbs = the_matrix_label == 'A' ? checkOtherMatrix('A', iRow, iCol, POrig) : checkOtherMatrix('P', iRow, iCol, AOrig);

    // -------------------------------------------------------------------------
    // EXCHANGE ACTION when initial useGibbs = true and check if Gibbs is usable.
    double s = 0, su = 0, mean = 0, sd = 0; // EJF -- MFO check
    std::pair<double, double> alphaparam(0,0);

    if (useGibbs)
    {
        alphaparam = GAPSNorm::calcAlphaParameters(the_matrix_label, _nFactor, D, S, AOrig, POrig, iGene1,
            iPattern1, iGene2, iPattern2, iSample1, iSample2);
        s = alphaparam.first * _annealingTemperature;
        su = alphaparam.second * _annealingTemperature;
        mean = su / s;
        sd = 1. / sqrt(s);
    }

    if (s == 0.0 && su == 0.0)
    {
        useGibbs = false;
    }

    // -------------------------------------------------------------------------
    if (useGibbs)
    {
        // set newMass1
        // need to retain exponential prior
        double plower = Random::p_norm(-mass1, mean, sd);
        double pupper = Random::p_norm(mass2, mean, sd);
        double u = plower + Random::uniform() * (pupper - plower);

        // must sample from prior if the computed parameters are not good for Gibbs
        if (!(plower >  0.95 || pupper < 0.05 || s < epsilon ||
        newMass1 == DOUBLE_POSINF || newMass1 == DOUBLE_NEGINF))
        {
            proposal.deltaMass1 = Random::q_norm(u, mean, sd);
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
    double pnew = mass1 > mass2 ? Random::d_gamma(newMass1, 2.0, 1.0 / lambda)
        : Random::d_gamma(newMass2, 2.0, 1.0 / lambda)
    
    double pold = newMass1 > newMass2 ? Random::d_gamma(mass1, 2.0, 1.0 / lambda)
        : Random::d_gamma(mass2, 2.0, 1.0 / lambda);

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

    if (1 - priorLL + delLLnew * _annealingTemperature  < log(Random::uniform()))
    {
        domain.acceptProposal(false);
        update_sysChi2(delLLnew);
        return true;
    }
    return false;
}

#endif