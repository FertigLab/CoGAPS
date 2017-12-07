#include "GibbsSampler.h"

#include "MatAlgo.h"

// ************* METHODS FOR COMPUTING LIKELIHOOD FUNCTIONS ******************
// ---------------------------------------------------------------------------
double GibbsSampler::cal_logLikelihood()
{
    return GAPSNorm::calChi2(_DMatrix, _SMatrix, _AMatrix, _PMatrix,
        _nFactor) / 2.0;
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

    unsigned int iRow = matrix[0].row;
    unsigned int iCol = matrix[0].col;

    double delLL = computeDeltaLL(the_matrix_label, D, S, AOrig, POrig, _matrixElemChange);
    domain.acceptProposal(proposal);
    matrix.elemUpdate(matChange);
    update_sysChi2(delLL);

    // rebirth
    double newMass = -proposal.deltaMass1;
    bool gibbs = the_matrix_label == 'A' ? checkOtherMatrix('A', iRow, iCol, POrig) : checkOtherMatrix('P', iRow, iCol, AOrig);

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

bool GibbsSampler::birth(char the_matrix_label, const Matrix &D, const Matrix &S,
const Matrix &AOrig, const Matrix &POrig, const AtomicProposal &proposal,
const AtomicSupport &domain)
{
    double origMass = proposal.deltaMass1;
    std::vector<ElementChange> matChange = domain.getElementChange(proposal);

    unsigned int iRow = matrix[0].row;
    unsigned int iCol = matrix[0].col;

    bool gibbs = the_matrix_label == 'A' ? checkOtherMatrix('A', iRow, iCol, POrig) : checkOtherMatrix('P', iRow, iCol, AOrig);

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

    // preparing quantities for possible Gibbs computation later.
    //EJF bool exchange = false;
    double priorLL = 0.;
    bool useGibbs = true;

    bool anyNonzero = false;
    switch (the_matrix_label)
    {
        case 'A':
            for (int jSample = 0; jSample < _nCol; jSample++)
            {
                anyNonzero = POrig(matChange[0].col,jSample) > epsilon || 
                    POrig(matChange[1].col,jSample) > epsilon;
                break;
            }
            break;

        case 'P':
            for (int jGene = 0; jGene < _nRow; jGene++)
            {
                anyNonzero = AOrig(jGene, matChange[0].row) > epsilon ||
                    AOrig(jGene, matChange[1],row) > epsilon;
                break;
            }
            break;
    }
    if (!anyNonzero)
    {
        useGibbs = false;
    }

    // -------------------------------------------------------------------------
    // EXCHANGE ACTION when initial useGibbs = true and check if Gibbs is usable.
    double s = 0, su = 0, mean = 0, sd = 0; // EJF -- MFO check
    std::pair<double, double> alphaparam(0,0);

    if (useGibbs == true)
    {
        if (domain.inDomain(loc1) && domain.inDomain(loc2))
        {
            alphaparam = GAPSNorm::calcAlphaParameters(the_matrix_label, _nFactor, D, S, AOrig, POrig, iGene1,
                iPattern1, iGene2, iPattern2, iSample1, iSample2);
            s = alphaparam.first;
            su = alphaparam.second;
            s = s * _annealingTemperature;
            su = su * _annealingTemperature;
            mean = su / s;
            sd = 1. / sqrt(s);
        }
    }

    if (s == 0.0 && su == 0.0)
    {
        useGibbs = false;
    }

    // -------------------------------------------------------------------------
    if (useGibbs == true)
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
            double gibbsMass1 = Random::q_norm(u, mean, sd);
            gibbsMass1 = std::max(gibbsMass1, -mass1);
            gibbsMass1 = std::min(gibbsMass1, mass2);

            double gibbsMass2 = -gibbsMass1;
            proposal.deltaMass1 = gibbsMass1;
            proposal.deltaMass2 = gibbsMass2;

            std::vector<ElementChange> matChange = domain.getElementChange(proposal);
            double delLLnew = computeDeltaLL(the_matrix_label, D, S, AOrig, POrig, matChange);
            domain.acceptProposal(false);
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

// ************ METHODS FOR LOOPING AND CONTROL *****************************
void GibbsSampler::set_iter(unsigned long ext_iter)
{
    _iter = ext_iter;
}

// -----------------------------------------------------------------------------
void GibbsSampler::set_AnnealingTemperature()
{
    double SASteps = _nEquil;
    double SATemp = std::min(((double)_iter + 1.0) / (SASteps / 2.0), 1.0);

    if (SATemp < 0.0)
    {
        throw std::logic_error("Invalid annealing temperature.");
    }
    _annealingTemperature = SATemp;
}

// -----------------------------------------------------------------------------
// TODO should be called internally when updated, not by the top level
void GibbsSampler::check_atomic_matrix_consistency(char label)
{
    double total_atom_mass = label == 'A' ? _AAtomicdomain.totalMass()
        : _PAtomicdomain.totalMass();
    double total_matrix_mass = label == 'A' ? MatAlgo::sum(_AMatrix)
        : MatAlgo::sum(_PMatrix);

    double diff_total_mass = fabs(total_atom_mass - total_matrix_mass);

    if (diff_total_mass > 0.00001)
    {
        throw std::logic_error("Mass inconsistency between atomic domain and matrix!");
    }
}

