#include <testthat.h>
#include "../testthat-tweak.h"
#include "../utils/GapsPrint.h"
#include "../GapsParameters.h"
#include "../math/Random.h"

#include "../gibbs_sampler/SingleThreadedGibbsSampler.h"
#include "../gibbs_sampler/DenseNormalModel.h"
#include "../gibbs_sampler/SparseNormalModel.h"

#define INIT_SAMPLER(name, Sampler, Data) Sampler<Data> name(initGibbsSampler<Sampler, Data>())

// create dummy data, introduce roughly 25% sparsity
Matrix getDummyData(unsigned nrow, unsigned ncol)
{
    Matrix data(nrow, ncol);
    for (unsigned i = 0; i < data.nRow(); ++i)
    {
        for (unsigned j = 0; j < data.nCol(); ++j)
        {
            data(i,j) = i * j % 2 == 1 ? 0.f : static_cast<float>(i * j);
        }
    }
    return data;
}

template <template<class> class GibbsSampler, class DataModel>
GibbsSampler<DataModel> initGibbsSampler()
{
    // initialization parameters
    Matrix data(getDummyData(25, 50));
    GapsParameters params(data);
    // GapsRng stores a pointer to GapsRandomState (const GapsRandomState *mRandState),
    // so randState must outlive the sampler returned by value. Make it static: one
    // instance per template instantiation, alive for the program's lifetime.
    // (data/params are copied into the model, only randState is held by pointer.)
    static GapsRandomState randState(params.seed);
    return GibbsSampler<DataModel>(data, false, false, 0.01f, 100.f, params, &randState);
}

TEST_CASE("Sampler Construction", "[samplerhighlevel][construction]")
{
    SECTION("Build all sampler variants with default uncertainty")
    {
    // construct samplers using default uncertainty
    INIT_SAMPLER(sampler1, SingleThreadedGibbsSampler, DenseNormalModel);
    INIT_SAMPLER(sampler2, SingleThreadedGibbsSampler, SparseNormalModel);

    REQUIRE(sampler1.dataSparsity() == sampler2.dataSparsity());
    } // closes SECTION
}

// construct an A/P sampler pair for the given DataModel, run sampling, and
// require the fit (chiSq) to improve -- exercises the full update() path
template <class DataModel>
static void requireUpdateDecreasesChiSq()
{
    Matrix data(getDummyData(25, 50));
    GapsRandomState randState(42);
    GapsParameters params(data);
    SingleThreadedGibbsSampler<DataModel> A(data, true, false, params.alphaA,
        params.maxGibbsMassA, params, &randState);
    SingleThreadedGibbsSampler<DataModel> P(data, false, false, params.alphaP,
        params.maxGibbsMassP, params, &randState);

    // chiSq() uses the other matrix, which is only valid after sync()
    A.sync(P); P.sync(A);
    A.extraInitialization(); P.extraInitialization();

    double AChiInit = A.chiSq();
    double PChiInit = P.chiSq();

    // interleave update+sync as the algorithm does (sparse rebuilds its lookup
    // tables in sync(), so each sampler must sync after the other matrix changes)
    A.update(100); P.sync(A);
    P.update(100); A.sync(P);
    A.extraInitialization(); P.extraInitialization();

    REQUIRE(A.chiSq() < AChiInit);
    REQUIRE(P.chiSq() < PChiInit);
}

TEST_CASE("Sampler Update", "[samplerhighlevel][update]")
{
    SECTION("update() improves the fit -- dense model")
    {
        requireUpdateDecreasesChiSq<DenseNormalModel>();
    }
    SECTION("update() improves the fit -- sparse model")
    {
        requireUpdateDecreasesChiSq<SparseNormalModel>();
    }
}
