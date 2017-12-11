#include "catch.h"
#include "../Matrix.h"
#include "../Algorithms.h"


TEST_CASE("Test Algorithms.h")
{
    unsigned nrow = 250;
    unsigned ncol = 500;

    Vector v(nrow);
    TwoWayMatrix D(nrow, ncol), S(nrow, ncol), AP(nrow, ncol);
    RowMatrix P(nrow, ncol);
    ColMatrix A(nrow, ncol);

    for (unsigned r = 0; r < nrow; ++r)
    {
        v(r) = r / 100.0;
        for (unsigned c = 0; c < ncol; ++c)
        {
            D.set(r, c, (r + c) / 100.0);
            S.set(r, c, (r + c) / 100.0);
            AP.set(r, c, (r - c) / 100.0);
            P(r,c) = (r * c) / 100.0;
            A(r,c) = (r * c) / 100.0;
        }
    }

    SECTION("Simple Algorithms")
    {
        REQUIRE(gaps::algo::sum(v) == 311.25);
        REQUIRE(gaps::algo::sum(D) == 467500);
        
        REQUIRE(gaps::algo::mean(D) == 3.74);

        REQUIRE(gaps::algo::nonZeroMean(D) > 3.74);
        REQUIRE(gaps::algo::nonZeroMean(D) < 3.75);

        REQUIRE(gaps::algo::isRowZero(P, 0));
        REQUIRE(!gaps::algo::isRowZero(P, 1));
        
        REQUIRE(gaps::algo::isColZero(A, 0));
        REQUIRE(!gaps::algo::isColZero(A, 1));
    }

    SECTION("delLL calculation")
    {
        REQUIRE(gaps::algo::deltaLL(MatrixChange('A', 3, 2, 2.5), D, S, A, P, AP)
            == Approx(-83583201873.509).epsilon(0.001));
    }

    SECTION("alpha parameters calculation")
    {
        AlphaParameters ap = gaps::algo::alphaParameters(
            MatrixChange('A', 3, 2, 2.5), D, S, A, P, AP);
        REQUIRE(ap.s == Approx(1887.023).epsilon(0.001));
        REQUIRE(ap.su == Approx(-33433278390.625).epsilon(0.001));
    }
}

