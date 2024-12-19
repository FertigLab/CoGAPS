#include <testthat.h>
#include "../testthat-tweak.h"
#include "../data_structures/SparseMatrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/TsvParser.h"
#include "../file_parser/MtxParser.h"
#include "../math/Random.h"
#include "../math/VectorMath.h"
#include "../math/MatrixMath.h"

static std::vector<unsigned> sequentialVector(unsigned n)
{
    std::vector<unsigned> vec;
    for (unsigned i = 1; i <= n; ++i) // mimic R indices
    {
        vec.push_back(i);
    }
    return vec;
}

template <class DataType>
static void testFullConstructor(float expectedSum, unsigned nr, unsigned nc, 
const DataType &data, bool transpose=false, bool partitionRows=false,
const std::vector<unsigned> &indices=std::vector<unsigned>())
{
    SparseMatrix mat(data, transpose, partitionRows, indices);

    REQUIRE(mat.nRow() == nr);
    REQUIRE(mat.nCol() == nc);
    REQUIRE(expectedSum == gaps::sum(mat));
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

TEST_CASE("Test Writing/Reading Sparse Matrices from File")
{
    // matrix to use for testing
    Matrix ref(25, 50);
    GapsRandomState randState(123);
    GapsRng rng(&randState);
    for (unsigned i = 0; i < ref.nRow(); ++i)
    {
        for (unsigned j = 0; j < ref.nCol(); ++j)
        {
            ref(i,j) = (i + j) * (rng.uniform() < 0.5f ? 0.f : 1.f);
        }
    }

    // write matrix to file
    FileParser::writeToTsv("testMatWrite.tsv", ref);
    FileParser::writeToCsv("testMatWrite.csv", ref);
    FileParser::writeToMtx("testMatWrite.mtx", ref);

    // read matrices from file
    SparseMatrix mat(ref, false, false, sequentialVector(0));
    SparseMatrix matTsv("testMatWrite.tsv", false, false, sequentialVector(0));
    SparseMatrix matCsv("testMatWrite.csv", false, false, sequentialVector(0));
    SparseMatrix matMtx("testMatWrite.mtx", false, false, sequentialVector(0));

    // delete files
    std::remove("testMatWrite.tsv");
    std::remove("testMatWrite.csv");
    std::remove("testMatWrite.mtx");

    // test matrices
    REQUIRE(gaps::sum(mat) == gaps::sum(ref));
    REQUIRE(gaps::sum(matTsv) == gaps::sum(ref));
    REQUIRE(gaps::sum(matCsv) == gaps::sum(ref));
    REQUIRE(gaps::sum(matMtx) == gaps::sum(ref));
}

TEST_CASE("Test SparseMatrix.h")
{
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
}
