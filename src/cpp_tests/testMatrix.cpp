#include <testthat.h>
#include "../Matrix.h"

#include <algorithm>
#include <cmath>

CATCH_TEST_CASE("Test Matrix.h")
{
    std::vector< std::vector<double> > data;
    data.push_back({ 1.0,  2.0,  3.0});
    data.push_back({-4.0,  5.0, -6.0});
    data.push_back({ 7.0, -8.0,  9.0});
    
    Matrix testMat (data, 'A');
    CATCH_REQUIRE(testMat.cal_mean() == 1.0);
    CATCH_REQUIRE(testMat.cal_totalsum() == 9.0);
}

