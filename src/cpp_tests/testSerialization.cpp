#include "catch.h"
#include "../utils/Archive.h"
#include "../data_structures/Matrix.h"
#include "../math/Random.h"
#include "../atomic/AtomicDomain.h"
#include "../atomic/ProposalQueue.h"

// put Archive in it's own scope so it gets destructed (file stream closed)

TEST_CASE("Reading/Writing to an Archive")
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
}

TEST_CASE("Serialization of primitive types")
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
}

TEST_CASE("Vector Serialization")
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
}

TEST_CASE("HybridVector Serialization")
{

}

TEST_CASE("SparseVector Serialization")
{

}

TEST_CASE("Matrix Serialization")
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
}

TEST_CASE("HybridMatrix Serialization")
{

}

TEST_CASE("SparseMatrix Serialization")
{

}

TEST_CASE("Random Generator Serialization")
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
}

TEST_CASE("GibbsSampler Serialization")
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

TEST_CASE("GapsParameters Serialization")
{

}

TEST_CASE("GapsStatistics Serialization")
{

}

#ifdef GAPS_INTERNAL_TESTS
TEST_CASE("AtomicDomain Serialization")
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

TEST_CASE("ProposalQueue Serialization")
{
    const unsigned nGenes = 10000;
    const unsigned nPatterns = 100;
    const unsigned nIterations = 1000;

    GapsRandomState randStateWrite(123);
    AtomicDomain domainWrite(nGenes * nPatterns);
    ProposalQueue queueWrite(nGenes, nPatterns, &randStateWrite);
    queueWrite.setAlpha(0.01f);
    queueWrite.setLambda(0.01f);

    for (unsigned i = 0; i < nIterations; ++i)
    {
        queueWrite.populate(domainWrite, nIterations);
        for (unsigned j = 0; j < queueWrite.size(); ++j)
        {
            switch (queueWrite[j].type)
            {
                case 'B':
                    queueWrite.acceptBirth();
                    queueWrite[j].atom1->mass = 3.f;
                    break;
                case 'D':
                    queueWrite.acceptDeath();
                    domainWrite.erase(queueWrite[j].atom1->pos);
                    break;
                case 'M':
                    queueWrite[j].atom1->pos = queueWrite[j].pos;
                    break;
                case 'E':
                    float mass1 = queueWrite[j].atom1->mass;
                    queueWrite[j].atom1->mass = queueWrite[j].atom2->mass;
                    queueWrite[j].atom2->mass = mass1;
                    break;
            }
        }
        queueWrite.clear();
    }

    {
        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << randStateWrite << queueWrite << domainWrite;
    }

    GapsRandomState randStateRead(456);
    AtomicDomain domainRead(nGenes * nPatterns);
    ProposalQueue queueRead(nGenes, nPatterns, &randStateRead);
    queueRead.setAlpha(100.f);
    queueRead.setLambda(100.f);
    
    {
        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> randStateRead >> queueRead >> domainRead;
    }

    GapsRng rngWriteTest(&randStateWrite);
    GapsRng rngReadTest(&randStateRead);

    REQUIRE(rngWriteTest.uniform() == rngReadTest.uniform());
    REQUIRE(rngWriteTest.uniform() == rngReadTest.uniform());
    REQUIRE(domainWrite.size() == domainRead.size());
    REQUIRE(queueWrite.size() == 0);
    REQUIRE(queueRead.size() == 0);

#ifdef GAPS_INTERNAL_TESTS
    REQUIRE(queueWrite.mRng.uniform() == queueRead.mRng.uniform());
    REQUIRE(queueWrite.mRng.uniform() == queueRead.mRng.uniform());
    REQUIRE(queueWrite.mRng.uniform() == queueRead.mRng.uniform());
    REQUIRE(queueWrite.mRng.uniform() == queueRead.mRng.uniform());
    REQUIRE(queueWrite.mMinAtoms == queueRead.mMinAtoms);
    REQUIRE(queueWrite.mMaxAtoms == queueRead.mMaxAtoms);
    REQUIRE(queueWrite.mBinLength == queueRead.mBinLength);
    REQUIRE(queueWrite.mNumCols == queueRead.mNumCols);
    REQUIRE(queueWrite.mAlpha == queueRead.mAlpha);
    REQUIRE(queueWrite.mDomainLength == queueRead.mDomainLength);
    REQUIRE(queueWrite.mNumBins == queueRead.mNumBins);
    REQUIRE(queueWrite.mUseCachedRng == queueRead.mUseCachedRng);
    for (unsigned i = 0; i < domainWrite.size(); ++i)
    {
        REQUIRE(domainWrite.mAtoms[i]->pos == domainRead.mAtoms[i]->pos);
        REQUIRE(domainWrite.mAtoms[i]->mass == domainRead.mAtoms[i]->mass);
    }
#endif

    for (unsigned i = 0; i < nIterations; ++i)
    {
        queueWrite.populate(domainWrite, nIterations);
        queueRead.populate(domainRead, nIterations);
        REQUIRE(queueWrite.size() == queueRead.size());

        for (unsigned j = 0; j < queueWrite.size(); ++j)
        {
            REQUIRE(queueWrite[j].type == queueRead[j].type);
            REQUIRE(queueWrite[j].pos == queueRead[j].pos);
            REQUIRE(queueWrite[j].r1 == queueRead[j].r1);
            REQUIRE(queueWrite[j].c1 == queueRead[j].c1);
            REQUIRE(queueWrite[j].r2 == queueRead[j].r2);
            REQUIRE(queueWrite[j].c2 == queueRead[j].c2);

            if (queueWrite[j].atom1 == NULL)
            {
                bool b = queueRead[j].atom1 == NULL;
                REQUIRE(b); // prevent pointer comparison warning in catch.h
            }
            else
            {
                REQUIRE(queueWrite[j].atom1->pos == queueRead[j].atom1->pos);
                REQUIRE(queueWrite[j].atom1->mass == queueRead[j].atom1->mass);
            }

            if (queueWrite[j].atom2 == NULL)
            {
                bool b = queueRead[j].atom2 == NULL;
                REQUIRE(b); // prevent pointer comparison warning in catch.h
            }
            else
            {
                REQUIRE(queueWrite[j].atom2->pos == queueRead[j].atom2->pos);
                REQUIRE(queueWrite[j].atom2->mass == queueRead[j].atom2->mass);
            }

            // process proposal
            switch (queueWrite[j].type)
            {
                case 'B':
                    queueWrite.acceptBirth();
                    queueRead.acceptBirth();
                    queueWrite[j].atom1->mass = 3.f;
                    queueRead[j].atom1->mass = 3.f;
                    break;
                case 'D':
                    queueWrite.acceptDeath();
                    queueRead.acceptDeath();
                    domainWrite.erase(queueWrite[j].atom1->pos);
                    domainRead.erase(queueRead[j].atom1->pos);
                    break;
                case 'M':
                    queueWrite[j].atom1->pos = queueWrite[j].pos;
                    queueRead[j].atom1->pos = queueRead[j].pos;
                    break;
                case 'E':
                    float mass1w = queueWrite[j].atom1->mass;
                    queueWrite[j].atom1->mass = queueWrite[j].atom2->mass;
                    queueWrite[j].atom2->mass = mass1w;

                    float mass1r = queueRead[j].atom1->mass;
                    queueRead[j].atom1->mass = queueRead[j].atom2->mass;
                    queueRead[j].atom2->mass = mass1r;
                    break;
            }
        }

        queueWrite.clear();
        queueRead.clear();
    }

    // cleanup directory
    std::remove("test_ar.temp");
}


