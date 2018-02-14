#ifndef __GAPS_GAPS_RUNNER_H__
#define __GAPS_GAPS_RUNNER_H__

// holds the data and dispatches the top level jobs
class GapsRunner
{
private:

    // Amplitude and Pattern matrices
    ColMatrix mAMatrix;
    RowMatrix mPMatrix;

    // used when updating A matrix
    RowMatrix mDMatrix;
    RowMatrix mSMatrix;
    RowMatrix mAPMatrix_A;

    // used when upating P matrix
    ColMatrix mDMatrix;
    ColMatrix mSMatrix;
    ColMatrix mAPMatrix_P;

    // gibbs sampler
    AmplitudeGibbsSampler mAGibbsSampler;
    PatternGibbsSampler mPGibbsSampler;

    // proposal queue
    ProposalQueue mAQueue;
    ProposalQueue mPQueue;
    
    // atomic domain
    AtomicDomain mADomain;
    AtomicDomina mPDomain;

    // number of cores available for jobs
    unsigned mNumCores;

public:

    GapsRunner() {}

    void run(unsigned nASteps, unsigned nPSteps)
    {
        update(mADomain, mAQueue, mAGibbsSampler, nASteps);
        mAPMatrix_P = mAPMatrix_A;

        update(mPDomain, mPQueue, mPGibbsSampler, nPSteps);
        mAPMatrix_A = mAPMatrix_P;
    }

    // Performance Metrics
    // 1) % of cores used in each iteration
    // 2) given fixed nCores, how does speed get better with matrix size,
    //      i.e. when does the overhead of parallelization start paying off
    // 3) % of program spent in parallel portion, i.e. not in populate queue
    void update(AtomicDomain domain, ProposalQueue queue, GibbsSampler sampler,
    unsigned nUpdates)
    {
        unsigned n = 0;
        while (n < nSteps)
        {
            // want this to be as quick as possible - otherwise there would be
            // a large speed up to making this run concurrently along with the
            // processProposal jobs, but that is much, much more complicated
            // to implement
            assert(nSteps - (queue.size() + n) >= 0);
            queue.populate(domain, nSteps - (queue.size() + n))

            unsigned nJobs = std::min(queue.size(), mNumCores);
            for (unsigned i = 0; i < nJobs; ++i) // can be run in parallel
            {
                sampler.processProposal(domain, queue[i]);
            }
            queue.clear(nJobs);
            n += nJobs;
            assert(n <= nSteps);
        }
    }
};

#endif