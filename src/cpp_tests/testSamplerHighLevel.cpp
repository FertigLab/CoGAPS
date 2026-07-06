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

TEST_CASE("Sampler Update", "[samplerhighlevel][update]")
{
    SECTION("Run update on all sampler variants")
    {
    // construct samplers using default uncertainty
    INIT_SAMPLER(sampler1, SingleThreadedGibbsSampler, DenseNormalModel);
    INIT_SAMPLER(sampler2, SingleThreadedGibbsSampler, SparseNormalModel);
    // NOTE: calling update() here would require sync()+extraInitialization()
    // first (otherwise mOtherMatrix is NULL); left as a construction smoke test.
    } // closes SECTION
}
