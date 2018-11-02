#include "catch.h"
#include "../gibbs_sampler/DenseGibbsSampler.h"
#include "../gibbs_sampler/SparseGibbsSampler.h"

#define TEST_APPROX(x) Approx(x).epsilon(0.001f)

TEST_CASE("Test SparseGibbsSampler")
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
        SparseGibbsSampler ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        SparseGibbsSampler PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);
    
        ASampler.sync(PSampler);
        PSampler.sync(ASampler);

        REQUIRE(ASampler.chiSq() == 100.f * data.nRow() * data.nCol());
        REQUIRE(PSampler.chiSq() == 100.f * data.nRow() * data.nCol());

    #ifdef GAPS_DEBUG
        REQUIRE(ASampler.internallyConsistent());
        REQUIRE(PSampler.internallyConsistent());
    #endif
    }

#ifdef GAPS_INTERNAL_TESTS
    SECTION("Test consistency between alpha parameters calculations")
    {
        // create the "data"
        Matrix data(100, 75);
        GapsRandomState randState(123);
        GapsRng rng(&randState);
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned j = 0; j < data.nCol(); ++j)
            {
                data(i,j) = rng.uniform32(1,14) * (rng.uniform() < 0.5f ? 0.f : 1.f);
            }
        }

        // create pair of sparse gibbs samplers
        GapsParameters params(data);
        SparseGibbsSampler sparse_ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        SparseGibbsSampler sparse_PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);
        sparse_ASampler.sync(sparse_PSampler);
        sparse_PSampler.sync(sparse_ASampler);

        // create pair of dense gibbs samplers
        DenseGibbsSampler dense_ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        DenseGibbsSampler dense_PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);
        dense_ASampler.sync(dense_PSampler);
        dense_PSampler.sync(dense_ASampler);

        // set the A and P matrix to the same thing
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned k = 0; k < params.nPatterns; ++k)
            {
                float val = rng.uniform(0.f, 10.f) * (rng.uniform() < 0.2f ? 0.f : 1.f);
                dense_ASampler.mMatrix(i,k) = val;
                sparse_ASampler.mMatrix.add(i, k, val);
            }
        }        
        REQUIRE(gaps::sum(dense_ASampler.mMatrix) == gaps::sum(sparse_ASampler.mMatrix));

        for (unsigned j = 0; j < data.nCol(); ++j)
        {
            for (unsigned k = 0; k < params.nPatterns; ++k)
            {
                float val = rng.uniform(0.f, 10.f) * (rng.uniform() < 0.2f ? 0.f : 1.f);
                dense_PSampler.mMatrix(j,k) = val;
                sparse_PSampler.mMatrix.add(j, k, val);
            }
        }        
        REQUIRE(gaps::sum(dense_PSampler.mMatrix) == gaps::sum(sparse_PSampler.mMatrix));

        // sync them back up
        sparse_ASampler.sync(sparse_PSampler);
        sparse_PSampler.sync(sparse_ASampler);
        dense_ASampler.sync(dense_PSampler);
        dense_PSampler.sync(dense_ASampler);
        dense_ASampler.extraInitialization();
        dense_PSampler.extraInitialization();

///////////////// test that alphaParameters are the same ///////////////////////
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned k = 0; k < params.nPatterns; ++k)
            {
                AlphaParameters sa = sparse_ASampler.alphaParameters(i,k);
                AlphaParameters da = dense_ASampler.alphaParameters(i,k);
                REQUIRE(sa.s >= 0.f);
                REQUIRE(da.s >= 0.f);
                if (sa.s <= gaps::epsilon || da.s <= gaps::epsilon)
                {
                    REQUIRE(sa.s <= gaps::epsilon);
                    REQUIRE(da.s <= gaps::epsilon);
                }
                REQUIRE(sa.s == TEST_APPROX(da.s));
                REQUIRE(sa.s_mu == TEST_APPROX(da.s_mu));
            }
        }

        for (unsigned j = 0; j < data.nCol(); ++j)
        {
            for (unsigned k = 0; k < params.nPatterns; ++k)
            {
                AlphaParameters sa = sparse_PSampler.alphaParameters(j,k);
                AlphaParameters da = dense_PSampler.alphaParameters(j,k);
                REQUIRE(sa.s >= 0.f);
                REQUIRE(da.s >= 0.f);
                if (sa.s <= gaps::epsilon || da.s <= gaps::epsilon)
                {
                    REQUIRE(sa.s <= gaps::epsilon);
                    REQUIRE(da.s <= gaps::epsilon);
                }
                REQUIRE(sa.s == TEST_APPROX(da.s));
                REQUIRE(sa.s_mu == TEST_APPROX(da.s_mu));
            }
        }

///////////// test two dimensional alphaParameters are the same ////////////////
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned k1 = 0; k1 < params.nPatterns; ++k1)
            {
                for (unsigned k2 = k1+1; k2 < params.nPatterns; ++k2)
                {
                    AlphaParameters sa = sparse_ASampler.alphaParameters(i,k1,i,k2);
                    AlphaParameters da = dense_ASampler.alphaParameters(i,k1,i,k2);
                    REQUIRE(sa.s >= 0.f);
                    REQUIRE(da.s >= 0.f);
                    if (sa.s <= gaps::epsilon || da.s <= gaps::epsilon)
                    {
                        REQUIRE(sa.s <= gaps::epsilon);
                        REQUIRE(da.s <= gaps::epsilon);
                    }
                    REQUIRE(sa.s == TEST_APPROX(da.s));
                    REQUIRE(sa.s_mu == TEST_APPROX(da.s_mu));

                    // symmetry
                    sa = sparse_ASampler.alphaParameters(i,k2,i,k1);
                    da = dense_ASampler.alphaParameters(i,k2,i,k1);
                    REQUIRE(sa.s >= 0.f);
                    REQUIRE(da.s >= 0.f);
                    if (sa.s <= gaps::epsilon || da.s <= gaps::epsilon)
                    {
                        REQUIRE(sa.s <= gaps::epsilon);
                        REQUIRE(da.s <= gaps::epsilon);
                    }
                    REQUIRE(sa.s == TEST_APPROX(da.s));
                    REQUIRE(sa.s_mu == TEST_APPROX(da.s_mu));
                }
            }
        }

        for (unsigned j = 0; j < data.nCol(); ++j)
        {
            for (unsigned k1 = 0; k1 < params.nPatterns; ++k1)
            {
                for (unsigned k2 = k1+1; k2 < params.nPatterns; ++k2)
                {
                    AlphaParameters sa = sparse_PSampler.alphaParameters(j,k1,j,k2);
                    AlphaParameters da = dense_PSampler.alphaParameters(j,k1,j,k2);
                    REQUIRE(sa.s >= 0.f);
                    REQUIRE(da.s >= 0.f);
                    if (sa.s <= gaps::epsilon || da.s <= gaps::epsilon)
                    {
                        REQUIRE(sa.s <= gaps::epsilon);
                        REQUIRE(da.s <= gaps::epsilon);
                    }
                    REQUIRE(sa.s == TEST_APPROX(da.s));
                    REQUIRE(sa.s_mu == TEST_APPROX(da.s_mu));

                    // symmetry
                    sa = sparse_PSampler.alphaParameters(j,k2,j,k1);
                    da = dense_PSampler.alphaParameters(j,k2,j,k1);
                    REQUIRE(sa.s >= 0.f);
                    REQUIRE(da.s >= 0.f);
                    if (sa.s <= gaps::epsilon || da.s <= gaps::epsilon)
                    {
                        REQUIRE(sa.s <= gaps::epsilon);
                        REQUIRE(da.s <= gaps::epsilon);
                    }
                    REQUIRE(sa.s == TEST_APPROX(da.s));
                    REQUIRE(sa.s_mu == TEST_APPROX(da.s_mu));
                }
            }
        }

///////////// test alphaParameters with change are the same ////////////////////
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned k = 0; k < params.nPatterns; ++k)
            {
                float ch = rng.uniform(0.f, 25.f);
                AlphaParameters sa = sparse_ASampler.alphaParametersWithChange(i,k,ch);
                AlphaParameters da = dense_ASampler.alphaParametersWithChange(i,k,ch);
                REQUIRE(sa.s >= 0.f);
                REQUIRE(da.s >= 0.f);
                if (sa.s <= gaps::epsilon || da.s <= gaps::epsilon)
                {
                    REQUIRE(sa.s <= gaps::epsilon);
                    REQUIRE(da.s <= gaps::epsilon);
                }
                REQUIRE(sa.s == TEST_APPROX(da.s));
                REQUIRE(sa.s_mu == TEST_APPROX(da.s_mu));
            }
        }

        for (unsigned j = 0; j < data.nCol(); ++j)
        {
            for (unsigned k = 0; k < params.nPatterns; ++k)
            {
                float ch = rng.uniform(0.f, 25.f);
                AlphaParameters sa = sparse_PSampler.alphaParametersWithChange(j,k,ch);
                AlphaParameters da = dense_PSampler.alphaParametersWithChange(j,k,ch);
                REQUIRE(sa.s >= 0.f);
                REQUIRE(da.s >= 0.f);
                if (sa.s <= gaps::epsilon || da.s <= gaps::epsilon)
                {
                    REQUIRE(sa.s <= gaps::epsilon);
                    REQUIRE(da.s <= gaps::epsilon);
                }
                REQUIRE(sa.s == TEST_APPROX(da.s));
                REQUIRE(sa.s_mu == TEST_APPROX(da.s_mu));
            }
        }
    }
#endif
}