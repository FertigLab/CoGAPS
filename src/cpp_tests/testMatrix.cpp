#include "catch.h"
#include "../data_structures/Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/TsvParser.h"
#include "../file_parser/MtxParser.h"

#if 0

TEST_CASE("Test Matrix.h")
{
    SECTION("Matrix/Vector Initialization")
    {
        Vector v(10);
        RowMatrix rm(10, 25);
        ColMatrix cm(10, 25);

        REQUIRE(v.size() == 10);
        REQUIRE(rm.nRow() == 10);
        REQUIRE(rm.nCol() == 25);
        REQUIRE(cm.nRow() == 10);
        REQUIRE(cm.nCol() == 25);
    }
}

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