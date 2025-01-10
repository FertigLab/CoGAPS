#include <testthat.h>
#include "../testthat-tweak.h"
#include "../gibbs_sampler/SingleThreadedGibbsSampler.h"
#include "../gibbs_sampler/DenseNormalModel.h"
#include "../data_structures/Matrix.h"

#include <Rcpp.h>

// convert R to C++ data type
// it is copied from CoGAPS.cpp not to change the maun headers
static Matrix convertRMatrix(const Rcpp::NumericMatrix &rmat)
{
    Matrix mat(rmat.nrow(), rmat.ncol());
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            mat(i,j) = rmat(i,j);
        }
    }
    return mat;
}


TEST_CASE("Test DenseGibbsSampler on random matrix","[densesinglesampler][randommat]")
{
    SECTION("Construct from random data matrix")
    {
        
        Matrix data(25, 50);
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned j = 0; j < data.nCol(); ++j)
            {
                data(i,j) = i + j + 1.f;
            }
        }

        GapsRandomState randState(42);
        GapsParameters params(data);
        SingleThreadedGibbsSampler<DenseNormalModel> ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        SingleThreadedGibbsSampler<DenseNormalModel> PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);
        

        REQUIRE(ASampler.chiSq() == 100.f * data.nRow() * data.nCol());
        REQUIRE(PSampler.chiSq() == 100.f * data.nRow() * data.nCol());
    
        double AChiInit=ASampler.chiSq();
        double PChiInit=PSampler.chiSq();

        ASampler.sync(PSampler);
        PSampler.sync(ASampler);
        ASampler.extraInitialization();
        PSampler.extraInitialization();

        REQUIRE(ASampler.chiSq() == AChiInit);
        REQUIRE(PSampler.chiSq() == PChiInit);

    #ifdef GAPS_DEBUG
        REQUIRE(ASampler.internallyConsistent());
        REQUIRE(PSampler.internallyConsistent());
    #endif
        ASampler.sync(PSampler);
        PSampler.sync(ASampler);
        ASampler.extraInitialization();
        PSampler.extraInitialization();

        REQUIRE(ASampler.chiSq() == AChiInit);
        REQUIRE(PSampler.chiSq() == PChiInit);
    }
}

TEST_CASE("Test DenseGibbsSampler on gist matrix","[densesinglesampler][gistmat]")
{

    SECTION("Construct from gist matrix")
    {
        Rcpp::Environment env = Rcpp::Environment::global_env();
        Rcpp::Function load("data");
        //R function is data()
        load("GIST");
        Matrix data=convertRMatrix(env["GIST.matrix"]);

        GapsRandomState randState(42);
        GapsParameters params(data);
        SingleThreadedGibbsSampler<DenseNormalModel> ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        SingleThreadedGibbsSampler<DenseNormalModel> PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);
    
        double AChiInit=ASampler.chiSq();
        double PChiInit=PSampler.chiSq();
    
        ASampler.sync(PSampler);
        PSampler.sync(ASampler);
        ASampler.extraInitialization();
        PSampler.extraInitialization();

        REQUIRE(ASampler.chiSq() == AChiInit);
        REQUIRE(PSampler.chiSq() == PChiInit);

    #ifdef GAPS_DEBUG
        REQUIRE(ASampler.internallyConsistent());
        REQUIRE(PSampler.internallyConsistent());
    #endif
    }
}
