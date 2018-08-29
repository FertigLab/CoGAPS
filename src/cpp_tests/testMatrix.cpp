#include "catch.h"
#include "../data_structures/Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/TsvParser.h"
#include "../file_parser/MtxParser.h"
#include "../math/Algorithms.h"

#if 0

static std::vector<unsigned> sequentialVector(unsigned n)
{
    std::vector<unsigned> vec;
    for (unsigned i = 0; i < n; ++i)
    {
        vec.push_back(i);
    }
    return vec;
}

template <class DataType>
static void testFullConstructor(float expectedSum, unsigned nr, unsigned nc, 
const DataType &data, bool transpose=false, bool partitionRows=false,
const std::vector<unsigned> &indices=std::vector<unsigned>(1))
{
    RowMatrix rm(data, transpose, partitionRows, indices);
    ColMatrix cm(data, transpose, partitionRows, indices);

    REQUIRE(rm.nRow() == nr);
    REQUIRE(rm.nCol() == nc);
    REQUIRE(cm.nRow() == nr);
    REQUIRE(cm.nCol() == nc);

    REQUIRE(expectedSum == gaps::algo::sum(rm));
    REQUIRE(expectedSum == gaps::algo::sum(cm));
}

template <class DataType>
static void testAllConstructorSituations(const DataType &data, unsigned nr,
unsigned nc, unsigned nIndices, float sum1, float sum2, float sum3)
{
    // No Transpose, No Subset
    testFullConstructor(sum1, nr, nc, data, false);

    // Transpose, No Subset
    testFullConstructor(sum1, nc, nr, data, true);

    // No Transpose, Subset Rows
    testFullConstructor(sum2, nIndices, nc, data, false, true,
        sequentialVector(nIndices));

    // Transpose, Subset Rows
    testFullConstructor(sum3, nIndices, nr, data, true, true,
        sequentialVector(nIndices));

    // No Transpose, Subset Columns
    testFullConstructor(sum3, nr, nIndices, data, false, false,
        sequentialVector(nIndices));

    // Transpose, Subset Columns
    testFullConstructor(sum2, nc, nIndices, data, true, false,
        sequentialVector(nIndices));
}

TEST_CASE("Test Writing/Reading Matrices from File")
{
    // matrix to use for testing
    Matrix ref(25, 50);
    for (unsigned i = 0; i < ref.nRow(); ++i)
    {
        for (unsigned j = 0; j < ref.nCol(); ++j)
        {
            ref(i,j) = i + j;
        }
    }

    // write matrix to file
    FileParser::writeToTsv("testMatWrite.tsv", ref);
    FileParser::writeToCsv("testMatWrite.csv", ref);
    FileParser::writeToMtx("testMatWrite.mtx", ref);

    // read matrices from file
    RowMatrix rmTsv("testMatWrite.tsv", false, false, sequentialVector(1));
    RowMatrix rmCsv("testMatWrite.csv", false, false, sequentialVector(1));
    RowMatrix rmMtx("testMatWrite.mtx", false, false, sequentialVector(1));
    ColMatrix cmTsv("testMatWrite.tsv", false, false, sequentialVector(1));
    ColMatrix cmCsv("testMatWrite.csv", false, false, sequentialVector(1));
    ColMatrix cmMtx("testMatWrite.mtx", false, false, sequentialVector(1));

    // delete files
    std::remove("testMatWrite.tsv");
    std::remove("testMatWrite.csv");
    std::remove("testMatWrite.mtx");

    // test matrices
    REQUIRE(gaps::algo::sum(rmTsv) == gaps::algo::sum(ref));
    REQUIRE(gaps::algo::sum(rmCsv) == gaps::algo::sum(ref));
    REQUIRE(gaps::algo::sum(rmMtx) == gaps::algo::sum(ref));
    REQUIRE(gaps::algo::sum(cmTsv) == gaps::algo::sum(ref));
    REQUIRE(gaps::algo::sum(cmCsv) == gaps::algo::sum(ref));
    REQUIRE(gaps::algo::sum(cmMtx) == gaps::algo::sum(ref));
}

TEST_CASE("Test Matrix.h")
{
    SECTION("Default Construction")
    {
        RowMatrix rm(10, 25);
        ColMatrix cm(10, 25);

        REQUIRE(rm.nRow() == 10);
        REQUIRE(rm.nCol() == 25);
        REQUIRE(cm.nRow() == 10);
        REQUIRE(cm.nCol() == 25);
    }

    SECTION("Copy Construction")
    {
        RowMatrix rm1(10, 25);
        ColMatrix cm1(rm1);

        REQUIRE(cm1.nRow() == 10);
        REQUIRE(cm1.nCol() == 25);

        RowMatrix rm2(10, 25);
        ColMatrix cm2(rm1);

        REQUIRE(rm2.nRow() == 10);
        REQUIRE(rm2.nCol() == 25);
    }

    SECTION("Full Constructor")
    {
        Matrix ref(10, 25);
        for (unsigned i = 0; i < ref.nRow(); ++i)
        {
            for (unsigned j = 0; j < ref.nCol(); ++j)
            {
                ref(i,j) = i + j;
            }
        }

        // write matrix to file
        FileParser::writeToTsv("testMatWrite.tsv", ref);
        FileParser::writeToCsv("testMatWrite.csv", ref);
        FileParser::writeToMtx("testMatWrite.mtx", ref);

        // test
        testAllConstructorSituations(ref, 10, 25, 5, 4125.f, 1750.f, 325.f);
        testAllConstructorSituations("testMatWrite.tsv", 10, 25, 5, 4125.f, 1750.f, 325.f);
        testAllConstructorSituations("testMatWrite.csv", 10, 25, 5, 4125.f, 1750.f, 325.f);
        testAllConstructorSituations("testMatWrite.mtx", 10, 25, 5, 4125.f, 1750.f, 325.f);

        // delete files
        std::remove("testMatWrite.tsv");
        std::remove("testMatWrite.csv");
        std::remove("testMatWrite.mtx");
    }

    SECTION("Arithmetic")
    {

    }
}


#endif