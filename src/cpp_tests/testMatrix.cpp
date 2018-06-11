#include "catch.h"
#include "../data_structures/Matrix.h"

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

    SECTION("Matrix Initialization from R Matrix")
    {
        Rcpp::Function asMatrix("as.matrix");
        Rcpp::Environment pkgEnv;
        pkgEnv = Rcpp::Environment::namespace_env("CoGAPS");
        Rcpp::NumericMatrix rD = asMatrix(pkgEnv.find("GIST.D"));
        Rcpp::NumericMatrix rS = asMatrix(pkgEnv.find("GIST.D"));

        REQUIRE(rD.nrow() == 1363);
        REQUIRE(rD.ncol() == 9);

        REQUIRE(rS.nrow() == 1363);
        REQUIRE(rS.ncol() == 9);

        RowMatrix rmD(rD);
        RowMatrix rmS(rS);

        ColMatrix cmD(rD);
        ColMatrix cmS(rS);

        REQUIRE(rmD.nRow() == 1363);
        REQUIRE(rmD.nCol() == 9);

        REQUIRE(rmS.nRow() == 1363);
        REQUIRE(rmS.nCol() == 9);

        REQUIRE(cmD.nRow() == 1363);
        REQUIRE(cmD.nCol() == 9);

        REQUIRE(cmS.nRow() == 1363);
        REQUIRE(cmS.nCol() == 9);
    }
}