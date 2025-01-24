#include <testthat.h>
#include "../testthat-tweak.h"
#include "../gibbs_sampler/SingleThreadedGibbsSampler.h"
#include "../gibbs_sampler/DenseNormalModel.h"
#include "../data_structures/Matrix.h"

#include <Rcpp.h>

// convert R to C++ data type
// it is copied from CoGAPS.cpp not to change the main headers
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


//copied from Matrix.cpp, changed: std::ostream instead of Archive
std::ostream& operator<<(std::ostream &ar, const Matrix &mat)
{
    ar << mat.nRow() <<" x "<< mat.nCol() << std::endl;
    for (unsigned i = 0; i < mat.nRow(); ++i) {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            ar << mat(i,j) << " ";
        }
        ar<<std::endl;
    }
    return ar;
}



TEST_CASE("Basic test on tiny matrix","[densesinglesampler][tinymat]")
{
    SECTION("Construct tiny matrix and do steps")
    {
        unsigned int patno=2;
        
        Matrix data(5,10);
        for (unsigned i = 0; i < data.nRow(); ++i)
        {
            for (unsigned j = 0; j < data.nCol(); ++j)
            {
                data(i,j) = i + j + 1.f; //nonzero
            }
        }

        GapsRandomState randState(42);
        GapsParameters params(data);
        params.nPatterns=patno;
        SingleThreadedGibbsSampler<DenseNormalModel> ASampler(data, true, false, params.alphaA,
            params.maxGibbsMassA, params, &randState);
        Matrix AMM(ASampler.MyMatrix());
        AMM.pad(5);
        ASampler.setMatrix(AMM);
        SingleThreadedGibbsSampler<DenseNormalModel> PSampler(data, false, false, params.alphaP,
            params.maxGibbsMassP, params, &randState);
        Matrix PMM(PSampler.MyMatrix());
        PMM.pad(2);
        PSampler.setMatrix(PMM);
        ASampler.sync(PSampler);
        PSampler.sync(ASampler);
        ASampler.extraInitialization();
        //actually, it is AP = A times P
        PSampler.extraInitialization();
        //actually, it is AP = A times P
        const Matrix & AAP=ASampler.APMatrix();
        const Matrix & PAP=PSampler.APMatrix();
        //just a ref
        std::cout<<std::fixed<<std::setprecision(0)<<"A:\n"<<ASampler.MyMatrix()<<"\nP\n"<<PSampler.MyMatrix()<<"\nA.AP\n"<<AAP<<"\nP.AP\n"<<PAP<<"\n";
        REQUIRE(gaps::sum(AAP) == 20 * data.nRow() * data.nCol());
        REQUIRE(gaps::sum(PAP) == 20 * data.nRow() * data.nCol());
    }
}



TEST_CASE("Test DenseGibbsSampler on random matrix","[densesinglesampler][randommat]")
{
    SECTION("Construct from random data matrix and do steps")
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
        float A_APSumInit=gaps::sum(ASampler.APMatrix());
        float P_APSumInit=gaps::sum(PSampler.APMatrix());
        float A_SumInit=gaps::sum(ASampler.MyMatrix());
        float P_SumInit=gaps::sum(PSampler.MyMatrix());
        std::cout<<std::fixed<<std::setprecision(3);
        for (unsigned i = 0; i < 2; ++i)
        {   
            std::cout<<"A: "<<" sum before="<<gaps::sum(ASampler.MyMatrix())<<" ";
            ASampler.update(1, 1);
            std::cout<<" sum after="<<gaps::sum(ASampler.MyMatrix())<<"\n";
            std::cout<<"P: "<<" sum before="<<gaps::sum(PSampler.MyMatrix())<<" ";
            PSampler.update(1, 1);
            std::cout<<" sum after="<<gaps::sum(PSampler.MyMatrix())<<"\n";

            ASampler.extraInitialization();
            PSampler.extraInitialization();
            
            std::cout<<" AP sums: "<<
                gaps::sum(ASampler.APMatrix())<<" and "<<
                gaps::sum(PSampler.APMatrix())<<std::endl;

            /*ASampler.sync(PSampler);
            PSampler.sync(ASampler);

            std::cout<<" AP sums after sync: "<<
                gaps::sum(ASampler.APMatrix())<<" and "<<
                gaps::sum(PSampler.APMatrix())<<std::endl;*/
            ASampler.extraInitialization();
            PSampler.extraInitialization();
            std::cout<<" AP sums after extra init: "<<
                gaps::sum(ASampler.APMatrix())<<" and "<<
                gaps::sum(PSampler.APMatrix())<<std::endl;
            
        }
        ASampler.sync(PSampler);
        PSampler.sync(ASampler);

        ASampler.extraInitialization();
        PSampler.extraInitialization();
        
        REQUIRE(ASampler.chiSq() < AChiInit);
        REQUIRE(PSampler.chiSq() < PChiInit);

        REQUIRE(gaps::sum(ASampler.APMatrix()) != A_APSumInit);
        REQUIRE(gaps::sum(PSampler.APMatrix()) != P_APSumInit);

        REQUIRE(gaps::sum(ASampler.MyMatrix()) != A_SumInit);
        REQUIRE(gaps::sum(PSampler.MyMatrix()) != P_SumInit);
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
        float A_APSumInit=gaps::sum(ASampler.APMatrix());
        float P_APSumInit=gaps::sum(PSampler.APMatrix());
        
        float A_SumInit=gaps::sum(ASampler.MyMatrix());
        float P_SumInit=gaps::sum(PSampler.MyMatrix());

        ASampler.update(100, 1);
        PSampler.update(100, 1);

        ASampler.sync(PSampler);
        PSampler.sync(ASampler);
        
        ASampler.extraInitialization();
        PSampler.extraInitialization();

        REQUIRE(ASampler.chiSq() < AChiInit);
        REQUIRE(PSampler.chiSq() < PChiInit);

        REQUIRE(gaps::sum(ASampler.APMatrix()) != A_APSumInit);
        REQUIRE(gaps::sum(PSampler.APMatrix()) != P_APSumInit);

        REQUIRE(gaps::sum(ASampler.MyMatrix()) != A_SumInit);
        REQUIRE(gaps::sum(PSampler.MyMatrix()) != P_SumInit);


    }
}
