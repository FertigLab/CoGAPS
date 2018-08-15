#include "catch.h"
#include "../Archive.h"
#include "../data_structures/Matrix.h"
#include "../GibbsSampler.h"
#include "../math/Random.h"

#if 0

TEST_CASE("Test Archive.h")
{
    SECTION("Reading/Writing to an Archive")
    {
        Archive ar1("test_ar.temp", ARCHIVE_WRITE);
        ar1 << 3;
        ar1.close();

        Archive ar2("test_ar.temp", ARCHIVE_READ);
        unsigned i = 0;
        ar2 >> i;
        REQUIRE(i == 3);
        ar2.close();
    }

    SECTION("Serialization of primitive types")
    {
        // test values
        unsigned u_read = 0, u_write = 456;
        uint32_t u32_read = 0, u32_write = 512;
        uint64_t u64_read = 0, u64_write = 0xAABBCCDDEE;
        float f_read = 0.f, f_write = 0.123542f;
        double d_read = 0., d_write = 0.54362;
        bool b_read = false, b_write = true;

        // write to archive
        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << u_write;
        arWrite << u32_write;
        arWrite << u64_write;
        arWrite << f_write;
        arWrite << d_write;
        arWrite << b_write;
        arWrite.close();

        // read from archive
        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> u_read;
        arRead >> u32_read;
        arRead >> u64_read;
        arRead >> f_read;
        arRead >> d_read;
        arRead >> b_read;
        arRead.close();

        // test that values are the same
        REQUIRE(u_read == u_write);
        REQUIRE(u32_read == u32_write);
        REQUIRE(u64_read == u64_write);
        REQUIRE(f_read == f_write);
        REQUIRE(d_read == d_write);
        REQUIRE(b_read == b_write);
    }
    
    SECTION("Vector Serialization")
    {
        Vector vec_read(100), vec_write(100);
        for (unsigned i = 0; i < 100; ++i)
        {
            vec_write[i] = gaps::random::uniform(0.0, 2.0);
        }

        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << vec_write;
        arWrite.close();

        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> vec_read;
        arRead.close();

        REQUIRE(vec_read.size() == vec_write.size());

        for (unsigned i = 0; i < 100; ++i)
        {
            REQUIRE(vec_read[i] == vec_write[i]);
        }
    }

    SECTION("Matrix Serialization")
    {
        RowMatrix rMat_read(100,100), rMat_write(100,100);
        ColMatrix cMat_read(100,100), cMat_write(100,100);

        for (unsigned i = 0; i < 100; ++i)
        {
            for (unsigned j = 0; j < 100; ++j)
            {
                rMat_write(i,j) = gaps::random::uniform(0.0, 2.0);
                cMat_write(i,j) = gaps::random::uniform(0.0, 2.0);
            }
        }

        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        arWrite << rMat_write;
        arWrite << cMat_write;
        arWrite.close();

        Archive arRead("test_ar.temp", ARCHIVE_READ);
        arRead >> rMat_read;
        arRead >> cMat_read;
        arRead.close();

        REQUIRE(rMat_read.nRow() == rMat_write.nRow());
        REQUIRE(rMat_read.nCol() == rMat_write.nCol());
        REQUIRE(cMat_read.nRow() == cMat_write.nRow());
        REQUIRE(cMat_read.nCol() == cMat_write.nCol());
    
        for (unsigned i = 0; i < 100; ++i)
        {
            for (unsigned j = 0; j < 100; ++j)
            {
                REQUIRE(rMat_read(i,j) == rMat_write(i,j));
                REQUIRE(cMat_read(i,j) == cMat_write(i,j));
            }
        }
    }

    SECTION("GibbsSampler Serialization")
    {
        //TODO
    }

    SECTION("Random Generator Serialization")
    {
        std::vector<float> randSequence;

        gaps::random::setSeed(123);
        volatile float burn_in = 0.0;
        for (unsigned i = 0; i < 1000; ++i)
        {
            burn_in = gaps::random::uniform(0,1);
        }
        REQUIRE(burn_in < 1);

        Archive arWrite("test_ar.temp", ARCHIVE_WRITE);
        gaps::random::save(arWrite);
        arWrite.close();

        for (unsigned i = 0; i < 1000; ++i)
        {
            randSequence.push_back(gaps::random::uniform());
            randSequence.push_back(gaps::random::uniform(0.3, 6.4));
            randSequence.push_back(gaps::random::exponential(5.5));
        }
    
        gaps::random::setSeed(11);
        Archive arRead("test_ar.temp", ARCHIVE_READ);
        gaps::random::load(arRead);
        arRead.close();

        for (unsigned i = 0; i < 1000; ++i)
        {
            REQUIRE(gaps::random::uniform() == randSequence[i++]);
            REQUIRE(gaps::random::uniform(0.3, 6.4) == randSequence[i++]);
            REQUIRE(gaps::random::exponential(5.5) == randSequence[i]);
        }
    }

    // cleanup directory
    std::remove("test_ar.temp");
}

#endif