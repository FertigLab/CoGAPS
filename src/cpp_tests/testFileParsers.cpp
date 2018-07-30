#include "catch.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/TsvParser.h"
#include "../file_parser/MtxParser.h"

#include "../data_structures/Matrix.h"

#include <Rcpp.h>

TEST_CASE("Test Parsers")
{
    Rcpp::Environment env = Rcpp::Environment::global_env();
    std::string csvPath = Rcpp::as<std::string>(env["gistCsvPath"]);
    std::string tsvPath = Rcpp::as<std::string>(env["gistTsvPath"]);
    std::string mtxPath = Rcpp::as<std::string>(env["gistMtxPath"]);

    SECTION("Test CsvParser")
    {
        CsvParser p(csvPath);
        REQUIRE(p.nRow() == 1363);
        REQUIRE(p.nCol() == 9);

        unsigned row = 0;
        unsigned col = 0;
        unsigned count = 0;
        while (p.hasNext())
        {
            MatrixElement e(p.getNext());
            REQUIRE(e.row == row);
            REQUIRE(e.col == col);

            ++count;
            ++col;
            if (col == 9) {
                ++row;
                col = 0;
            }
        }
        REQUIRE(count == 12267);
    }

    SECTION("Test TsvParser")
    {
        TsvParser p(tsvPath);
        REQUIRE(p.nRow() == 1363);
        REQUIRE(p.nCol() == 9);

        unsigned row = 0;
        unsigned col = 0;
        unsigned count = 0;
        while (p.hasNext())
        {
            MatrixElement e(p.getNext());
            REQUIRE(e.row == row);
            REQUIRE(e.col == col);

            ++count;
            ++col;
            if (col == 9) {
                ++row;
                col = 0;
            }
        }
        REQUIRE(count == 12267);
    }

    SECTION("Test MtxParser")
    {
        MtxParser p(mtxPath);
        REQUIRE(p.nRow() == 1363);
        REQUIRE(p.nCol() == 9);

        unsigned count = 0;
        while (p.hasNext())
        {
            MatrixElement e(p.getNext());

            REQUIRE(e.row < 1363);
            REQUIRE(e.col < 9);
            ++count;
        }
        REQUIRE(count == 12267);
    }
}
