class GapsRunner
{
private:

    // atomic domain
    // proposal queue
    // matrices
    // statistics

public:

    // restrict inner loop to number of cores, then re-populate before calling again
    //  simple solution that might be hard to beat
    void update(unsigned nSteps)
    {
        queue.populate();
        unsigned size = std::min(queue.size(), nCores);
        for (unsigned i = 0; i < size; ++i)
        {
            processProposal(queue[i]); // dispatch to correct function
        }
    }

};

void processProposal(AtomicProposal proposal)
{
    switch (proposal.type)
    {
        case 'D': death(proposal.pos1); break;
        case 'B': birth(proposal.pos1); break;
        case 'M': move(proposal.pos1, proposal.pos2); break;
        case 'E': exchange(proposal.pos1, proposal.pos2); break;
    }
}

// A matrix

void birth(uint64_t pos, float mass)
{
    unsigned row = getRow(pos);
    unsigned col = getCol(pos);
    bool useGibbs = !gaps::algo::isRowZero(mPMatrix, row);
    float mass = useGibbs ? gibbsMass(row, col) : gaps::random::exponential(mLambda);
    domain.addAtom(pos, mass);
    matrix(row,col) += mass;
    updateAPMatrix(row, col, mass);
}

void death(uint64_t pos, float mass)
{
    unsigned row = getRow(pos);
    unsigned col = getCol(pos);
    matrix(row,col) += -mass;
    updateAPMatrix(row, col, -mass);

    bool useGibbs = !gaps::algo::isRowZero(mPMatrix, row);

    float newMass = useGibbs ? gibbsMass(row, col) : mass;
    float deltaLL = computeDeltaLL(row, col, newMass);
    if (deltaLL * temp >= std::log(gaps::random::uniform()))
    {
        matrix(row,col) += newMass;
        updateAPMatrix(row, col, newMass);
        domain.updateAtomMass(pos, newMass - mass);
        queue.minAtoms++;
    }
    else
    {
        domain.deleteAtom(pos);
        queue.maxAtoms--;
    }
}

// automatically accept moves/exchanges in the same bin

void move(uint64_t src, uint64_t dest)
{
    unsigned srcRow = getRow(src);
    unsigned srcCol = getCol(src);
    float mass = domain.at(src);
    unsigned destRow = getRow(dest);
    unsigned destCol = getCol(dest);
    float deltaLL = computeDeltaLL(srcRow, srcCol, -mass, destRow, destCol, mass);
    if (deltaLL * temp >= std::log(gaps::random::uniform()))
    {
        matrix(srcRow, srcCol) += -mass;
        matrix(destRow, destCol) += mass;
        updateAPMatrix(srcRow, srcCol, -mass);
        updateAPMatrix(destRow, destCol, mass);
        domain.deleteAtom(src);
        domain.addAtom(dest, mass);
    }
}

void exchange(uint64_t p1, uint64_t p2)
{
    if (!gaps::algo::isRowZero(mPMatrix, getRow(p1)) && !gaps::algo::isRowZero(mPmatrix, getRow(p2)))
    {
        AlphaParameters alphaParam = gaps::algo::alphaParameters();
        alphaParam.s *= mAnnealingTemp;
        alphaParam.su *= mAnnealingTemp;

        if (alphaParam.s > EPSILON)
        {
            float mean = alphaParam.su / alphaParam.s;
            float sd = 1.f / std::sqrt(alphaParam.s);
            float plower = gaps::random::p_norm(-mass1, mean, sd);
            float pupper = gaps::random::p_norm(mass2, mean, sd);

            if (!(plower >  0.95f || pupper < 0.05f))
            {
                float u = gaps::random::uniform(plower, pupper);
                prop.delta1 = gaps::random::q_norm(u, mean, sd);
        
                float delta1 = gaps::random::q_norm(u, mean, sd);
                delta1 = std::max(-mass1, std::min(delta1, mass2));
                float delta2 = -delta1;
            }
        }
    }

    float mass1 = domain.at(p1);
    float mass2 = domain.at(p2);
    float newMass1 = mass1 + delta1;
    float newMass2 = mass2 + delta2;

    float pnewMass = mass1 > mass2 ? newMass1 : newMass2;
    float poldMass = newMass1 > newMass2 ? mass1 : mass2;

    float pnew = gaps::random::d_gamma(pnewMass, 2.0, 1.f / domain.lambda());
    float pold = gaps::random::d_gamma(poldMass, 2.0, 1.f / domain.lambda());

    if (pold == 0.f && pnew != 0.f)
    {
        evaluateChange(domain, prop, change, 0.f, true);
    }
    else
    {
        float priorLL = (pold == 0.f) ? 0.f : log(pnew / pold);
        float rejectProb = std::log(gaps::random::uniform()) - priorLL;
        evaluateChange(domain, prop, change, rejectProb);
    }
    delta1 = domain.updateAtomMass(p1, delta1);
    delta2 = domain.updateAtomMass(p2, delta2);

    matrix(getRow(p1), getCol(p1)) += delta1;
    matrix(getRow(p2), getCol(p2)) += delta2;
    updateAPMatrix(getRow(p1), getCol(p1), delta1);
    updateAPMatrix(getRow(p2), getCol(p2), delta2);

    queue.maxAtoms = (domain.at(p1) == 0) ? queue.maxAtoms - 1 : queue.maxAtoms;
    queue.minAtoms = (domain.at(p1) == 0) ? queue.minAtoms : queue.minAtoms + 1;
    queue.maxAtoms = (domain.at(p2) == 0) ? queue.maxAtoms - 1 : queue.maxAtoms;
    queue.minAtoms = (domain.at(p2) == 0) ? queue.minAtoms : queue.minAtoms + 1;
}