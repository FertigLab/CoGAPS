#include "catch.h"
#include "../gibbs_sampler/GibbsSampler.h"
#include "../gibbs_sampler/DenseStoragePolicy.h"

TEST_CASE("Test DenseGibbsSampler")
{
    SECTION("Construct from data matrix")
    {
        Matrix data(25, 50);
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned j = 0; j < data.nCol(); ++j)
            {
                data(i,j) = i + j + 1.f;
            }
        }

        GapsRandomState randState(123);
        GapsParameters params(data);
        GibbsSampler<DenseStorage> ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        GibbsSampler<DenseStorage> PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);
    
        REQUIRE(ASampler.chiSq() == 100.f * data.nRow() * data.nCol());
        REQUIRE(PSampler.chiSq() == 100.f * data.nRow() * data.nCol());
    
        ASampler.sync(PSampler);
        PSampler.sync(ASampler);
        ASampler.extraInitialization();
        PSampler.extraInitialization();

        REQUIRE(ASampler.chiSq() == 100.f * data.nRow() * data.nCol());
        REQUIRE(PSampler.chiSq() == 100.f * data.nRow() * data.nCol());

    #ifdef GAPS_DEBUG
        REQUIRE(ASampler.internallyConsistent());
        REQUIRE(PSampler.internallyConsistent());
    #endif
    }
}