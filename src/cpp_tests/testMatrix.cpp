#include "catch.h"
#include "../data_structures/Matrix.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/TsvParser.h"
#include "../file_parser/MtxParser.h"

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

    SECTION("Matrix Initialization from .tsv file")
    {
        std::vector<unsigned> whichRows;
        std::vector<unsigned> whichCols;
        for (unsigned i = 0; i < 1363; ++i)
        {
            whichRows.push_back(i);
        }
        for (unsigned i = 0; i < 9; ++i)
        {
            whichCols.push_back(i);
        }
        TsvParser p1("../../inst/extdata/GIST.tsv");
        RowMatrix tsvRowMatrix(p1, true, whichRows);
        REQUIRE(tsvRowMatrix.nRow() == 1363);
        REQUIRE(tsvRowMatrix.nCol() == 9);

        TsvParser p2("../../inst/extdata/GIST.tsv");
        ColMatrix tsvColMatrix(p2, true, whichRows);
        REQUIRE(tsvColMatrix.nRow() == 1363);
        REQUIRE(tsvColMatrix.nCol() == 9);

        TsvParser p3("../../inst/extdata/GIST.tsv");
        RowMatrix tsvRowMatrixT(p3, false, whichCols);
        REQUIRE(tsvRowMatrixT.nRow() == 9);
        REQUIRE(tsvRowMatrixT.nCol() == 1363);

        TsvParser p4("../../inst/extdata/GIST.tsv");
        ColMatrix tsvColMatrixT(p4, false, whichCols);
        REQUIRE(tsvColMatrixT.nRow() == 9);
        REQUIRE(tsvColMatrixT.nCol() == 1363);

        float sum = 0;
        for (unsigned i = 0; i < tsvRowMatrix.nCol(); ++i)
        {
            sum += tsvRowMatrix.getRow(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < tsvColMatrix.nCol(); ++i)
        {
            sum += tsvColMatrix.getCol(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < tsvRowMatrixT.nRow(); ++i)
        {
            sum += tsvRowMatrixT.getRow(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < tsvRowMatrixT.nRow(); ++i)
        {
            sum += tsvColMatrixT.getCol(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        for (unsigned i = 0; i < 1000; ++i)
        {
            whichRows.pop_back();
        }

        for (unsigned i = 0; i < 7; ++i)
        {
            whichCols.pop_back();
        }

        TsvParser p5("../../inst/extdata/GIST.tsv");
        tsvRowMatrix = RowMatrix(p5, true, whichRows);
        REQUIRE(tsvRowMatrix.nRow() == 363);
        REQUIRE(tsvRowMatrix.nCol() == 9);

        TsvParser p6("../../inst/extdata/GIST.tsv");
        tsvColMatrix = ColMatrix(p6, true, whichRows);
        REQUIRE(tsvColMatrix.nRow() == 363);
        REQUIRE(tsvColMatrix.nCol() == 9);

        TsvParser p7("../../inst/extdata/GIST.tsv");
        tsvRowMatrixT = RowMatrix(p7, false, whichCols);
        REQUIRE(tsvRowMatrixT.nRow() == 2);
        REQUIRE(tsvRowMatrixT.nCol() == 1363);

        TsvParser p8("../../inst/extdata/GIST.tsv");
        tsvColMatrixT = ColMatrix(p8, false, whichCols);
        REQUIRE(tsvColMatrixT.nRow() == 2);
        REQUIRE(tsvColMatrixT.nCol() == 1363);

        sum = 0;
        for (unsigned i = 0; i < tsvRowMatrix.nCol(); ++i)
        {
            sum += tsvRowMatrix.getRow(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < tsvColMatrix.nCol(); ++i)
        {
            sum += tsvColMatrix.getCol(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < tsvRowMatrixT.nRow(); ++i)
        {
            sum += tsvRowMatrixT.getRow(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 20);

        sum = 0;
        for (unsigned i = 0; i < tsvRowMatrixT.nRow(); ++i)
        {
            sum += tsvColMatrixT.getCol(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 20);
    }

    SECTION("Matrix Initialization from .mtx file")
    {
        std::vector<unsigned> whichRows;
        std::vector<unsigned> whichCols;
        for (unsigned i = 0; i < 1363; ++i)
        {
            whichRows.push_back(i);
        }
        for (unsigned i = 0; i < 9; ++i)
        {
            whichCols.push_back(i);
        }

        MtxParser p1("../../inst/extdata/GIST.mtx");
        RowMatrix mtxRowMatrix(p1, true, whichRows);
        REQUIRE(mtxRowMatrix.nRow() == 1363);
        REQUIRE(mtxRowMatrix.nCol() == 9);

        MtxParser p2("../../inst/extdata/GIST.mtx");
        ColMatrix mtxColMatrix(p2, true, whichRows);
        REQUIRE(mtxColMatrix.nRow() == 1363);
        REQUIRE(mtxColMatrix.nCol() == 9);

        MtxParser p3("../../inst/extdata/GIST.mtx");
        RowMatrix mtxRowMatrixT = RowMatrix(p3, false, whichCols);
        REQUIRE(mtxRowMatrixT.nRow() == 9);
        REQUIRE(mtxRowMatrixT.nCol() == 1363);

        MtxParser p4("../../inst/extdata/GIST.mtx");
        ColMatrix mtxColMatrixT = ColMatrix(p4, false, whichCols);
        REQUIRE(mtxColMatrixT.nRow() == 9);
        REQUIRE(mtxColMatrixT.nCol() == 1363);

        float sum = 0;
        for (unsigned i = 0; i < mtxRowMatrix.nCol(); ++i)
        {
            sum += mtxRowMatrix.getRow(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < mtxColMatrix.nCol(); ++i)
        {
            sum += mtxColMatrix.getCol(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < mtxRowMatrixT.nRow(); ++i)
        {
            sum += mtxRowMatrixT.getRow(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < mtxRowMatrixT.nRow(); ++i)
        {
            sum += mtxColMatrixT.getCol(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        for (unsigned i = 0; i < 1000; ++i)
        {
            whichRows.pop_back();
        }

        for (unsigned i = 0; i < 7; ++i)
        {
            whichCols.pop_back();
        }

        MtxParser p5("../../inst/extdata/GIST.mtx");
        mtxRowMatrix = RowMatrix(p5, true, whichRows);
        REQUIRE(mtxRowMatrix.nRow() == 363);
        REQUIRE(mtxRowMatrix.nCol() == 9);

        MtxParser p6("../../inst/extdata/GIST.mtx");
        mtxColMatrix = ColMatrix(p6, true, whichRows);
        REQUIRE(mtxColMatrix.nRow() == 363);
        REQUIRE(mtxColMatrix.nCol() == 9);

        MtxParser p7("../../inst/extdata/GIST.mtx");
        mtxRowMatrixT = RowMatrix(p7, false, whichCols);
        REQUIRE(mtxRowMatrixT.nRow() == 2);
        REQUIRE(mtxRowMatrixT.nCol() == 1363);

        MtxParser p8("../../inst/extdata/GIST.mtx");
        mtxColMatrixT = ColMatrix(p8, false, whichCols);
        REQUIRE(mtxColMatrixT.nRow() == 2);
        REQUIRE(mtxColMatrixT.nCol() == 1363);

        sum = 0;
        for (unsigned i = 0; i < mtxRowMatrix.nCol(); ++i)
        {
            sum += mtxRowMatrix.getRow(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < mtxColMatrix.nCol(); ++i)
        {
            sum += mtxColMatrix.getCol(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < mtxRowMatrixT.nRow(); ++i)
        {
            sum += mtxRowMatrixT.getRow(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 20);

        sum = 0;
        for (unsigned i = 0; i < mtxRowMatrixT.nRow(); ++i)
        {
            sum += mtxColMatrixT.getCol(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 20);
    }

    SECTION("Matrix Initialization from .csv file")
    {
        std::vector<unsigned> whichRows;
        std::vector<unsigned> whichCols;
        for (unsigned i = 0; i < 1363; ++i)
        {
            whichRows.push_back(i);
        }
        for (unsigned i = 0; i < 9; ++i)
        {
            whichCols.push_back(i);
        }
        CsvParser p1("../../inst/extdata/GIST.csv");
        RowMatrix csvRowMatrix(p1, true, whichRows);
        REQUIRE(csvRowMatrix.nRow() == 1363);
        REQUIRE(csvRowMatrix.nCol() == 9);

        CsvParser p2("../../inst/extdata/GIST.csv");
        ColMatrix csvColMatrix(p2, true, whichRows);
        REQUIRE(csvColMatrix.nRow() == 1363);
        REQUIRE(csvColMatrix.nCol() == 9);

        CsvParser p3("../../inst/extdata/GIST.csv");
        RowMatrix csvRowMatrixT = RowMatrix(p3, false, whichCols);
        REQUIRE(csvRowMatrixT.nRow() == 9);
        REQUIRE(csvRowMatrixT.nCol() == 1363);

        CsvParser p4("../../inst/extdata/GIST.csv");
        ColMatrix csvColMatrixT = ColMatrix(p4, false, whichCols);
        REQUIRE(csvColMatrixT.nRow() == 9);
        REQUIRE(csvColMatrixT.nCol() == 1363);

        float sum = 0;
        for (unsigned i = 0; i < csvRowMatrix.nCol(); ++i)
        {
            sum += csvRowMatrix.getRow(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < csvColMatrix.nCol(); ++i)
        {
            sum += csvColMatrix.getCol(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < csvRowMatrixT.nRow(); ++i)
        {
            sum += csvRowMatrixT.getRow(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < csvRowMatrixT.nRow(); ++i)
        {
            sum += csvColMatrixT.getCol(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        for (unsigned i = 0; i < 1000; ++i)
        {
            whichRows.pop_back();
        }

        for (unsigned i = 0; i < 7; ++i)
        {
            whichCols.pop_back();
        }

        CsvParser p5("../../inst/extdata/GIST.csv");
        csvRowMatrix = RowMatrix(p5, true, whichRows);
        REQUIRE(csvRowMatrix.nRow() == 363);
        REQUIRE(csvRowMatrix.nCol() == 9);

        CsvParser p6("../../inst/extdata/GIST.csv");
        csvColMatrix = ColMatrix(p6, true, whichRows);
        REQUIRE(csvColMatrix.nRow() == 363);
        REQUIRE(csvColMatrix.nCol() == 9);

        CsvParser p7("../../inst/extdata/GIST.csv");
        csvRowMatrixT = RowMatrix(p7, false, whichCols);
        REQUIRE(csvRowMatrixT.nRow() == 2);
        REQUIRE(csvRowMatrixT.nCol() == 1363);

        CsvParser p8("../../inst/extdata/GIST.csv");
        csvColMatrixT = ColMatrix(p8, false, whichCols);
        REQUIRE(csvColMatrixT.nRow() == 2);
        REQUIRE(csvColMatrixT.nCol() == 1363);

        sum = 0;
        for (unsigned i = 0; i < csvRowMatrix.nCol(); ++i)
        {
            sum += csvRowMatrix.getRow(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < csvColMatrix.nCol(); ++i)
        {
            sum += csvColMatrix.getCol(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 97);

        sum = 0;
        for (unsigned i = 0; i < csvRowMatrixT.nRow(); ++i)
        {
            sum += csvRowMatrixT.getRow(i)[0];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 20);

        sum = 0;
        for (unsigned i = 0; i < csvRowMatrixT.nRow(); ++i)
        {
            sum += csvColMatrixT.getCol(0)[i];
        }
        sum *= 10;
        sum = (int) sum;
        REQUIRE(sum == 20);
    }
}
