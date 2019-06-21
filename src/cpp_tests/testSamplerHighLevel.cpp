#include "catch.h"
#include "../utils/GapsPrint.h"
#include "../GapsParameters.h"
#include "../math/Random.h"

#include "../gibbs_sampler/AsynchronousGibbsSampler.h"
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

template <class GibbsSampler>
template <class DataModel>
GibbsSampler<DataModel> initGibbsSampler()
{
    // initialization parameters
    Matrix data(getDummyData(25, 50));
    GapsParameters params(data);
    GapsRandomState randState(params.seed);
    return GibbsSampler<DataModel>(data, false, false, 0.01f, 100.f, params, randState);
}

TEST_CASE("Sampler Construction")
{
    // construct samplers using default uncertainty
    INIT_SAMPLER(sampler1, SingleThreadedGibbsSampler, DenseNormalModel);
    INIT_SAMPLER(sampler2, SingleThreadedGibbsSampler, SparseNormalModel);
    INIT_SAMPLER(sampler3, AsynchronousGibbsSampler, DenseNormalModel);
    INIT_SAMPLER(sampler4, AsynchronousGibbsSampler, SparseNormalModel);

    REQUIRE(sampler1.dataSparsity() == sampler2.dataSparsity());
    REQUIRE(sampler2.dataSparsity() == sampler3.dataSparsity());
    REQUIRE(sampler3.dataSparsity() == sampler4.dataSparsity());
}

TEST_CASE("Sampler Update")
{
    // construct samplers using default uncertainty
    INIT_SAMPLER(sampler1, SingleThreadedGibbsSampler, DenseNormalModel);
    INIT_SAMPLER(sampler2, SingleThreadedGibbsSampler, SparseNormalModel);
    INIT_SAMPLER(sampler3, AsynchronousGibbsSampler, DenseNormalModel);
    INIT_SAMPLER(sampler4, AsynchronousGibbsSampler, SparseNormalModel);
}