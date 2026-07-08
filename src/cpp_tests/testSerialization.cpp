#include <testthat.h>
#include "../testthat-tweak.h"
#include "../utils/Archive.h"
#include "../data_structures/Matrix.h"
#include "../data_structures/HybridVector.h"
#include "../data_structures/SparseVector.h"
#include "../data_structures/HybridMatrix.h"
#include "../data_structures/SparseMatrix.h"
#include "../math/Random.h"
#include "../math/MatrixMath.h"
#include "../atomic/AtomicDomain.h"
#include "../GapsParameters.h"
#include "../GapsStatistics.h"
#include "../gibbs_sampler/SingleThreadedGibbsSampler.h"
#include "../gibbs_sampler/DenseNormalModel.h"

// put Archive in it's own scope so it gets destructed (file stream closed)

TEST_CASE("Reading/Writing to an Archive", "[serialization][archive]")
{
    SECTION("Write an integer and read it back")
    {
    {
        Archive ar1("test_ar.temp", ARCHIVE_WRITE);
        ar1 << 3;
    } // in it's own scope so the file resources gets released

    {
        Archive ar2("test_ar.temp", ARCHIVE_READ);
        unsigned i = 0;
        ar2 >> i;
        REQUIRE(i == 3);
    }

    // cleanup directory
    std::remove("test_ar.temp");
    } // closes SECTION
}

TEST_CASE("Serialization of primitive types", "[serialization][primitives]")
{
    SECTION("Round-trip read/write of all primitive types")
    {
    // test values
    unsigned u_read = 0, u_write = 456;
    uint32_t u32_read = 0, u32_write = 512;
    uint64_t u64_read = 0, u64_write = 0xAABBCCDDEE;
    float f_read = 0.f, f_write = 0.123542f;
    double d_read = 0., d_write = 0.54362;
    bool b_read = false, b_write = true;

    // write to archive
    {
        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << u_write;
        arWrite << u32_write;
        arWrite << u64_write;
        arWrite << f_write;
        arWrite << d_write;
        arWrite << b_write;
    }

    // read from archive
    {
        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> u_read;
        arRead >> u32_read;
        arRead >> u64_read;
        arRead >> f_read;
        arRead >> d_read;
        arRead >> b_read;
    }

    // test that values are the same
    REQUIRE(u_read == u_write);
    REQUIRE(u32_read == u32_write);
    REQUIRE(u64_read == u64_write);
    REQUIRE(f_read == f_write);
    REQUIRE(d_read == d_write);
    REQUIRE(b_read == b_write);

    // cleanup directory
    std::remove("test_ar.temp");
    } // closes SECTION
}

TEST_CASE("Vector Serialization", "[serialization][vector]")
{
    SECTION("Round-trip read/write of Vector")
    {
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    Vector vec_read(100), vec_write(100);
    for (unsigned i = 0; i < 100; ++i)
    {
        vec_write[i] = rng.uniform(0.f, 2.f);
    }

    {
        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << vec_write;
    }


    {
        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> vec_read;
    }

    REQUIRE(vec_read.size() == vec_write.size());

    for (unsigned i = 0; i < 100; ++i)
    {
        REQUIRE(vec_read[i] == vec_write[i]);
    }

    // cleanup directory
    std::remove("test_ar.temp");
    } // closes SECTION
}

TEST_CASE("HybridVector Serialization", "[serialization][hybridvector]")
{
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    std::vector<float> in_v;
    for (unsigned n = 0; n < 100; ++n)
        in_v.push_back(n % 3 == 0 ? 0.f : rng.uniform(0.5f, 2.f));
    HybridVector vecWrite(in_v), vecRead(100);

    { Archive ar("test_ar.temp", ARCHIVE_WRITE); ar << vecWrite; }
    { Archive ar("test_ar.temp", ARCHIVE_READ);  ar >> vecRead; }

    REQUIRE(vecRead.size() == vecWrite.size());
    for (unsigned i = 0; i < vecWrite.size(); ++i)
        REQUIRE(vecRead[i] == vecWrite[i]);

    std::remove("test_ar.temp");
}

TEST_CASE("SparseVector Serialization", "[serialization][sparsevector]")
{
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    std::vector<float> in_v;
    for (unsigned n = 0; n < 100; ++n)
        in_v.push_back(n % 4 == 0 ? rng.uniform(0.5f, 2.f) : 0.f);
    SparseVector vecWrite(in_v), vecRead(100);

    { Archive ar("test_ar.temp", ARCHIVE_WRITE); ar << vecWrite; }
    { Archive ar("test_ar.temp", ARCHIVE_READ);  ar >> vecRead; }

    REQUIRE(vecRead.size() == vecWrite.size());
    Vector denseWrite(vecWrite.getDense()), denseRead(vecRead.getDense());
    for (unsigned i = 0; i < denseWrite.size(); ++i)
        REQUIRE(denseRead[i] == denseWrite[i]);

    std::remove("test_ar.temp");
}

TEST_CASE("Matrix Serialization", "[serialization][matrix]")
{
    SECTION("Round-trip read/write of Matrix")
    {
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    Matrix matRead(100,100), matWrite(100,100);

    for (unsigned i = 0; i < 100; ++i)
    {
        for (unsigned j = 0; j < 100; ++j)
        {
            matWrite(i,j) = rng.uniform(0.f, 2.f);
        }
    }

    {
        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << matWrite;
    }

    {
        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> matRead;
    }

    REQUIRE(matRead.nRow() == matWrite.nRow());
    REQUIRE(matRead.nCol() == matWrite.nCol());

    for (unsigned i = 0; i < 100; ++i)
    {
        for (unsigned j = 0; j < 100; ++j)
        {
            REQUIRE(matRead(i,j) == matWrite(i,j));
        }
    }

    // cleanup directory
    std::remove("test_ar.temp");
    } // closes SECTION
}

TEST_CASE("HybridMatrix Serialization", "[serialization][hybridmatrix]")
{
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    HybridMatrix matWrite(20, 12), matRead(20, 12);
    for (unsigned i = 0; i < 20; ++i)
        for (unsigned j = 0; j < 12; ++j)
            matWrite.set(i, j, (i + j) % 3 == 0 ? 0.f : rng.uniform(0.5f, 2.f));

    { Archive ar("test_ar.temp", ARCHIVE_WRITE); ar << matWrite; }
    { Archive ar("test_ar.temp", ARCHIVE_READ);  ar >> matRead; }

    REQUIRE(matRead.nRow() == matWrite.nRow());
    REQUIRE(matRead.nCol() == matWrite.nCol());
    for (unsigned i = 0; i < 20; ++i)
        for (unsigned j = 0; j < 12; ++j)
            REQUIRE(matRead(i,j) == matWrite(i,j));

    std::remove("test_ar.temp");
}

TEST_CASE("SparseMatrix Serialization", "[serialization][sparsematrix]")
{
    GapsRandomState randState(123);
    GapsRng rng(&randState);

    Matrix dense(30, 10);
    for (unsigned i = 0; i < 30; ++i)
        for (unsigned j = 0; j < 10; ++j)
            dense(i,j) = (i % 3 == 0) ? rng.uniform(0.5f, 2.f) : 0.f;
    SparseMatrix matWrite(dense, false, false, std::vector<unsigned>());
    // build the read target from an all-zero matrix so its columns start empty:
    // the round-trip must repopulate them (catches the SparseVector count bug)
    Matrix zeros(30, 10);
    SparseMatrix matRead(zeros, false, false, std::vector<unsigned>());

    { Archive ar("test_ar.temp", ARCHIVE_WRITE); ar << matWrite; }
    { Archive ar("test_ar.temp", ARCHIVE_READ);  ar >> matRead; }

    REQUIRE(matRead.nRow() == matWrite.nRow());
    REQUIRE(matRead.nCol() == matWrite.nCol());
    for (unsigned j = 0; j < matWrite.nCol(); ++j)
    {
        Vector colWrite(matWrite.getCol(j).getDense());
        Vector colRead(matRead.getCol(j).getDense());
        REQUIRE(colRead.size() == colWrite.size());
        for (unsigned i = 0; i < colWrite.size(); ++i)
            REQUIRE(colRead[i] == colWrite[i]);
    }

    std::remove("test_ar.temp");
}

TEST_CASE("Random Generator Serialization", "[serialization][random]")
{
    SECTION("Round-trip read/write of GapsRandomState")
    {
    std::vector<float> randSequence;

    GapsRandomState randStateWrite(123);
    GapsRng rngWrite(&randStateWrite);

    volatile float burn_in = 0.0;
    for (unsigned i = 0; i < 1000; ++i)
    {
        burn_in = rngWrite.uniform(0,1);
    }
    REQUIRE(burn_in < 1);

    {
        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << rngWrite;
    }
    
    for (unsigned i = 0; i < 1000; ++i)
    {
        randSequence.push_back(rngWrite.uniform());
        randSequence.push_back(rngWrite.uniform(0.3, 6.4));
        randSequence.push_back(rngWrite.exponential(5.5));
    }

    GapsRandomState randStateRead(456);
    GapsRng rngRead(&randStateRead);

    {
        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> rngRead;
    }

    for (unsigned i = 0; i < 1000; ++i)
    {
        REQUIRE(rngRead.uniform() == randSequence[i++]);
        REQUIRE(rngRead.uniform(0.3, 6.4) == randSequence[i++]);
        REQUIRE(rngRead.exponential(5.5) == randSequence[i]);
    }

    // cleanup directory
    std::remove("test_ar.temp");
    } // closes SECTION
}

TEST_CASE("GibbsSampler Serialization", "[serialization][gibbssampler]")
{
#if 0
    Rcpp::Environment env = Rcpp::Environment::global_env();
    std::string csvPath = Rcpp::as<std::string>(env["gistCsvPath"]);

    GibbsSampler Asampler(csvPath, false, 7, false, std::vector<unsigned>());
    GibbsSampler Psampler(csvPath, true, 7, false, std::vector<unsigned>());
    Asampler.sync(Psampler);
    Psampler.sync(Asampler);
    
    Asampler.update(10000, 1);

    Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
    arWrite << Asampler;
    arWrite.close();

    GibbsSampler savedAsampler(csvPath, false, 7, false, std::vector<unsigned>());
    Archive arRead("test_ar.temp", ARCHIVE_READ);
    arRead >> savedAsampler;
    arRead.close();

    // cleanup directory
    std::remove("test_ar.temp");
#endif
}

TEST_CASE("GapsParameters Serialization", "[serialization][gapsparameters]")
{
    Matrix data(40, 15);
    GapsParameters pWrite(data);
    pWrite.seed = 777;
    pWrite.nGenes = 40;
    pWrite.nSamples = 15;
    pWrite.nPatterns = 5;
    pWrite.nIterations = 321;
    pWrite.alphaA = 0.02f;
    pWrite.alphaP = 0.03f;
    pWrite.maxGibbsMassA = 50.f;
    pWrite.maxGibbsMassP = 60.f;
    pWrite.useSparseOptimization = true;
    pWrite.checkpointInterval = 42;

    { Archive ar("test_ar.temp", ARCHIVE_WRITE); ar << pWrite; }
    GapsParameters pRead(data);
    { Archive ar("test_ar.temp", ARCHIVE_READ);  ar >> pRead; }

    // these are exactly the fields the serialization operators read/write; a
    // mismatch (or a << / >> field-order slip) fails here
    REQUIRE(pRead.seed == pWrite.seed);
    REQUIRE(pRead.nGenes == pWrite.nGenes);
    REQUIRE(pRead.nSamples == pWrite.nSamples);
    REQUIRE(pRead.nPatterns == pWrite.nPatterns);
    REQUIRE(pRead.nIterations == pWrite.nIterations);
    REQUIRE(pRead.alphaA == pWrite.alphaA);
    REQUIRE(pRead.alphaP == pWrite.alphaP);
    REQUIRE(pRead.maxGibbsMassA == pWrite.maxGibbsMassA);
    REQUIRE(pRead.maxGibbsMassP == pWrite.maxGibbsMassP);
    REQUIRE(pRead.useSparseOptimization == pWrite.useSparseOptimization);
    REQUIRE(pRead.checkpointInterval == pWrite.checkpointInterval);

    std::remove("test_ar.temp");
}

TEST_CASE("GapsStatistics Serialization", "[serialization][gapsstatistics]")
{
    // build two dense samplers with real state so the statistics matrices are
    // non-trivial, then round-trip the accumulated statistics
    Matrix data(12, 8);
    data.pad(15.f);
    GapsRandomState randState(42);
    GapsParameters params(data);
    params.nPatterns = 3;

    SingleThreadedGibbsSampler<DenseNormalModel> A(data, true, false, params.alphaA,
        params.maxGibbsMassA, params, &randState);
    SingleThreadedGibbsSampler<DenseNormalModel> P(data, false, false, params.alphaP,
        params.maxGibbsMassP, params, &randState);
    A.sync(P); P.sync(A);
    A.extraInitialization(); P.extraInitialization();
    A.update(300); P.update(300);
    A.sync(P); P.sync(A);

    GapsStatistics statsWrite(data.nRow(), data.nCol(), params.nPatterns);
    for (unsigned i = 0; i < 5; ++i) statsWrite.update(A, P);
    Matrix ameanW(statsWrite.Amean()), pmeanW(statsWrite.Pmean());
    REQUIRE(gaps::sum(ameanW) > 0.f); // there is real state to serialize

    { Archive ar("test_ar.temp", ARCHIVE_WRITE); ar << statsWrite; }
    GapsStatistics statsRead(data.nRow(), data.nCol(), params.nPatterns);
    { Archive ar("test_ar.temp", ARCHIVE_READ);  ar >> statsRead; }

    Matrix ameanR(statsRead.Amean()), pmeanR(statsRead.Pmean());
    REQUIRE(ameanR.nRow() == ameanW.nRow());
    REQUIRE(ameanR.nCol() == ameanW.nCol());
    for (unsigned i = 0; i < ameanW.nRow(); ++i)
        for (unsigned j = 0; j < ameanW.nCol(); ++j)
            REQUIRE(ameanR(i,j) == ameanW(i,j));
    for (unsigned i = 0; i < pmeanW.nRow(); ++i)
        for (unsigned j = 0; j < pmeanW.nCol(); ++j)
            REQUIRE(pmeanR(i,j) == pmeanW(i,j));

    std::remove("test_ar.temp");
}

#if 0
TEST_CASE("AtomicDomain Serialization", "[serialization][atomicdomain]")
{
    GapsRandomState randState(123);
    GapsRng rng(&randState);
    
    AtomicDomain domainWrite(100000);

    for (unsigned i = 0; i < 1000; ++i)
    {
        domainWrite.insert(rng.uniform64(), rng.uniform(0.f, 100.f));
    }

    {
        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << domainWrite;
    }

    AtomicDomain domainRead(1);
    
    {
        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> domainRead;
    }

    REQUIRE(domainWrite.front()->pos == domainRead.front()->pos);
    REQUIRE(domainWrite.front()->mass == domainRead.front()->mass);
    REQUIRE(domainWrite.size() == domainRead.size());
    REQUIRE(domainWrite.mDomainLength == domainRead.mDomainLength);

    for (unsigned i = 0; i < domainWrite.size(); ++i)
    {
        REQUIRE(domainWrite.mAtoms[i]->pos == domainRead.mAtoms[i]->pos);
        REQUIRE(domainWrite.mAtoms[i]->mass == domainRead.mAtoms[i]->mass);
    }

    // cleanup directory
    std::remove("test_ar.temp");
}
#endif

// Regression: SingleThreadedGibbsSampler's operator>> read the DataModel with >>
// but then chained << (write) for mDomain/mNumBins/.../mAlpha, so a checkpoint
// restore did NOT reload the atomic domain (a resumed run diverged from a fresh
// one). Fixed to use >>. This round-trip catches it: nAtoms() of the restored
// sampler must match the saved one.
TEST_CASE("GibbsSampler serialization round-trip","[serialization][gibbssampler-roundtrip]")
{
    Matrix data(12, 8);
    data.pad(15.f);
    GapsRandomState randState(42);
    GapsParameters params(data);
    params.nPatterns = 3;

    SingleThreadedGibbsSampler<DenseNormalModel> A(data, true, false, params.alphaA,
        params.maxGibbsMassA, params, &randState);
    SingleThreadedGibbsSampler<DenseNormalModel> P(data, false, false, params.alphaP,
        params.maxGibbsMassP, params, &randState);
    A.sync(P); P.sync(A);
    A.extraInitialization(); P.extraInitialization();
    A.update(500);
    REQUIRE(A.nAtoms() > 0); // there is domain state worth serializing

    {
        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << A;
    }

    SingleThreadedGibbsSampler<DenseNormalModel> Aread(data, true, false, params.alphaA,
        params.maxGibbsMassA, params, &randState);
    REQUIRE(Aread.nAtoms() == 0); // freshly constructed: empty domain
    {
        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> Aread;
    }

    // with the << bug the domain was never restored, so nAtoms() stayed 0
    REQUIRE(Aread.nAtoms() == A.nAtoms());
    REQUIRE(gaps::sum(Aread.MyMatrix()) == gaps::sum(A.MyMatrix()));

    std::remove("test_ar.temp");
}
