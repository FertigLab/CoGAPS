#include "catch.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/TsvParser.h"
#include "../file_parser/MtxParser.h"

#include "../data_structures/Matrix.h"

#include <Rcpp.h>

TEST_CASE("Test Parsers")
{
    SECTION("Test CsvParser")
    {
        CsvParser p("../../inst/extdata/GIST.csv");
        REQUIRE(p.nRow() == 1363);
        REQUIRE(p.nCol() == 9);

        unsigned row = 0;
        unsigned col = 0;
        unsigned count = 0;
        while (p.hasNext())
        {
            MatrixElement e(p.getNext());
            REQUIRE(e.row() == row);
            REQUIRE(e.col() == col);

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
        TsvParser p("../../inst/extdata/GIST.tsv");
        REQUIRE(p.nRow() == 1363);
        REQUIRE(p.nCol() == 9);

        unsigned row = 0;
        unsigned col = 0;
        unsigned count = 0;
        while (p.hasNext())
        {
            MatrixElement e(p.getNext());
            REQUIRE(e.row() == row);
            REQUIRE(e.col() == col);

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
        MtxParser p("../../inst/extdata/GIST.mtx");
        REQUIRE(p.nRow() == 1363);
        REQUIRE(p.nCol() == 9);

        unsigned row = 0;
        unsigned col = 0;
        unsigned count = 0;
        while (p.hasNext())
        {
            MatrixElement e(p.getNext());
            REQUIRE(e.row() == row);
            REQUIRE(e.col() == col);

            ++count;
            ++row;
            if (row == 1363) {
                ++col;
                row = 0;
            }
        }
        REQUIRE(count == 12267);
    }
}
