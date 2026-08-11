#include <testthat.h>
#include <Rcpp.h>
#include "../testthat-tweak.h"
#include "../file_parser/CharacterDelimitedParser.h"
#include "../file_parser/MtxParser.h"
#include "../data_structures/Matrix.h"

TEST_CASE("Test Parsers", "[fileparsers]")
{
    Rcpp::Function systemfile("system.file");
    std::string csvPath = Rcpp::as<std::string>(systemfile("extdata", "GIST.csv", Rcpp::Named("package") = "CoGAPS"));
    std::string tsvPath = Rcpp::as<std::string>(systemfile("extdata", "GIST.tsv", Rcpp::Named("package") = "CoGAPS"));
    std::string mtxPath = Rcpp::as<std::string>(systemfile("extdata", "GIST.mtx", Rcpp::Named("package") = "CoGAPS"));
    std::string gctPath = Rcpp::as<std::string>(systemfile("extdata", "GIST.gct", Rcpp::Named("package") = "CoGAPS"));

    SECTION("Test CsvParser")
    {
        CharacterDelimitedParser p(csvPath, ',');
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
            if (col == 9)
            {
                ++row;
                col = 0;
            }
        }
        REQUIRE(count == 12267);
    }

    SECTION("Test TsvParser")
    {
        CharacterDelimitedParser p(tsvPath, '\t');
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
            if (col == 9)
            {
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

    SECTION("Test GctParser")
    {
        CharacterDelimitedParser p(gctPath, '\t', true);
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
            if (col == 9)
            {
                ++row;
                col = 0;
            }
        }
        REQUIRE(count == 12267);
    }
}
