#include "catch.h"
#include "../data_structures/Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/TsvParser.h"
#include "../file_parser/MtxParser.h"
#include "../math/Algorithms.h"

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
static void testFullConstructor(unsigned nr, unsigned nc, float expectedSum,
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
static void testAllConstructorSituations(const DataType &data)
{
    // No Transpose, No Subset
    testFullConstructor(0.f, 10, 25, data, false);

    // Transpose, No Subset
    testFullConstructor(0.f, 25, 10, data, true);

    // No Transpose, Subset Rows
    testFullConstructor(0.f, 5, 25, data, false, true, sequentialVector(5))

    // Transpose, Subset Rows
    testFullConstructor(0.f, 5, 10, data, true, true, sequentialVector(5))

    // No Transpose, Subset Columns
    testFullConstructor(0.f, 10, 5, data, false, false, sequentialVector(5))

    // Transpose, Subset Columns
    testFullConstructor(0.f, 25, 5, data, true, false, sequentialVector(5))
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

    Matrix ref(10, 25);
    for (unsigned i = 0; i < ref.nRow(); ++i)
    {
        for (unsigned j = 0; j < ref.nCol(); ++j)
        {
            ref(i,j) = i + j;
        }
    }

    testAllConstructorSituations(ref);
    testAllConstructorSituations(ref);
    testAllConstructorSituations(ref);
    testAllConstructorSituations(ref);    


    SECTION("Construct from File - No Subset")
    {

    }

    SECTION("Construct from File - Subset")
    {

    }

    SECTION("Assignment")
    {

    }

    SECTION("Get Row/Col")
    {

    }

    SECTION("arithmetic")
    {

    }
}

#if 0

static void populateSequential(std::vector<unsigned> &vec, unsigned n)
{
    for (unsigned i = 0; i < n; ++i)
    {
        vec.push_back(i);
    }
}

static void testMatrixConstruction(const std::string &path)
{
    std::vector<unsigned> allRows, allCols;
    std::vector<unsigned> someRows, someCols;

    populateSequential(allRows, 1363);
    populateSequential(allCols, 9);
    populateSequential(someRows, 363);
    populateSequential(someCols, 2);

    FileParser p1(path);
    RowMatrix allRowMat(p1, true, allRows);
    REQUIRE(allRowMat.nRow() == 1363);
    REQUIRE(allRowMat.nCol() == 9);

    FileParser p2(path);
    ColMatrix allColMat(p2, true, allRows);
    REQUIRE(allColMat.nRow() == 1363);
    REQUIRE(allColMat.nCol() == 9);

    FileParser p3(path);
    RowMatrix allRowMatT(p3, false, allCols);
    REQUIRE(allRowMatT.nRow() == 9);
    REQUIRE(allRowMatT.nCol() == 1363);

    FileParser p4(path);
    ColMatrix allColMatT(p4, false, allCols);
    REQUIRE(allColMatT.nRow() == 9);
    REQUIRE(allColMatT.nCol() == 1363);

    for (unsigned i = 0; i < 1363; ++i)
    {
        REQUIRE(allRowMat(i, 3) == allColMat(i, 3));
        REQUIRE(allRowMat(i, 3) == allRowMatT(3, i));
        REQUIRE(allRowMat(i, 3) == allColMatT(3, i));
    }

/*
    FileParser p5(path);
    RowMatrix someRowMat(p1, true, someRows);
    REQUIRE(someRowMat.nRow() == 363);
    REQUIRE(someRowMat.nCol() == 9);

    FileParser p6(path);
    ColMatrix someColMat(p2, true, someRows);
    REQUIRE(someColMat.nRow() == 363);
    REQUIRE(someColMat.nCol() == 9);

    FileParser p7(path);
    RowMatrix someRowMatT(p3, false, someCols);
    REQUIRE(someRowMatT.nRow() == 2);
    REQUIRE(someRowMatT.nCol() == 1363);

    FileParser p8(path);
    ColMatrix someColMatT(p4, false, someCols);
    REQUIRE(someColMatT.nRow() == 2);
    REQUIRE(someColMatT.nCol() == 1363);

    for (unsigned i = 0; i < 363; ++i)
    {
        REQUIRE(someRowMat(i, 1) == someColMat(i, 1));
        REQUIRE(someRowMat(i, 1) == someRowMatT(1, i));
        REQUIRE(someRowMat(i, 1) == someColMatT(1, i));
    }
*/
}

#ifdef __GAPS_R_BUILD__
TEST_CASE("Test Matrix Construction from file")
{
    Rcpp::Environment env = Rcpp::Environment::global_env();
    std::string csvPath = Rcpp::as<std::string>(env["gistCsvPath"]);
    std::string tsvPath = Rcpp::as<std::string>(env["gistTsvPath"]);
    std::string mtxPath = Rcpp::as<std::string>(env["gistMtxPath"]);

    testMatrixConstruction(csvPath);
    testMatrixConstruction(tsvPath);
    testMatrixConstruction(mtxPath);
}
#endif

#endif