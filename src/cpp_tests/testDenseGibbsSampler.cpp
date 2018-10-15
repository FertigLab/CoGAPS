#include "catch.h"
#include "../gibbs_sampler/DenseGibbsSampler.h"

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

        GapsParameters params(data);
        DenseGibbsSampler ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params);
        DenseGibbsSampler PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params);
    
        REQUIRE(ASampler.chiSq() == 100.f * data.nRow() * data.nCol());
        REQUIRE(PSampler.chiSq() == 100.f * data.nRow() * data.nCol());
    
        ASampler.sync(PSampler);
        PSampler.sync(ASampler);
        ASampler.recalculateAPMatrix();
        PSampler.recalculateAPMatrix();

        REQUIRE(ASampler.chiSq() == 100.f * data.nRow() * data.nCol());
        REQUIRE(PSampler.chiSq() == 100.f * data.nRow() * data.nCol());

    #ifdef GAPS_DEBUG
        REQUIRE(ASampler.internallyConsistent());
        REQUIRE(PSampler.internallyConsistent());
    #endif
    }
}