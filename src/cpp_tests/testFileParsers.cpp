#include "catch.h"
#include "../file_parser/CsvParser.h"
#include "../file_parser/TsvParser.h"
#include "../file_parser/MtxParser.h"

#include "../data_structures/Matrix.h"

#include <Rcpp.h>

TEST_CASE("Test Parsers")
{
    SECTION("Test RowMatrix")
    {    
        //CsvParser csv("data/GIST.csv");
        //TsvParser tsv("data/GIST.tsv");
        //MtxParser mtx("data/GIST.mtx");

        //Rcpp::Environment pkgEnv;
        //pkgEnv = Rcpp::Environment::namespace_env("CoGAPS");
        //std::string mtxPath = pkgEnv.find("gistMtxPath");

        Rcpp::Function systemFile("system.file");
        //Rcpp::string mtxPath = systemFile("data/GIST.mtx", "CoGAPS");

        std::ifstream is("/mnt/c/Users/tsherma4/Documents/CoGAPS/Repo/data/GIST.mtx");
        std::string line;
        std::getline(is, line);
        std::cout << "\n" <<  line << "\nTHIS IS TEST OUTPUT\n";
        is.close();
    }
}