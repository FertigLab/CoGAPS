#include "catch.h"
#include "../Matrix.h"
#include "../Algorithms.h"


TEST_CASE("Test Algorithms.h")
{
    unsigned nrow = 25;
    unsigned ncol = 10;

    Vector v(nrow);
    TwoWayMatrix D(nrow, ncol), S(nrow, ncol), AP(nrow, ncol);
    RowMatrix P(nrow, ncol);
    ColMatrix A(nrow, ncol);

    for (unsigned r = 0; r < nrow; ++r)
    {
        v[r] = r;
        for (unsigned c = 0; c < ncol; ++c)
        {
            D.set(r, c, r + c);
            S.set(r, c, (r + c) / 100.f);
            AP.set(r, c, r - c);
            P(r,c) = r * c;
            A(r,c) = r * c;
        }
    }

    SECTION("sum")
    {
        REQUIRE(gaps::algo::sum(v) == 300);
        REQUIRE(gaps::algo::sum(D) == 300 * 10 + 45 * 25);
        REQUIRE(gaps::algo::sum(D) == gaps::algo::sum(S) * 100.0);
        REQUIRE(gaps::algo::sum(A) == gaps::algo::sum(P));
    }

    SECTION("mean")
    {
        float dTotal = 300 * 10 + 45 * 25;
    
        REQUIRE(gaps::algo::mean(D) == gaps::algo::mean(S) * 100.f);
        REQUIRE(gaps::algo::mean(D) == dTotal / (nrow * ncol));
        REQUIRE(gaps::algo::nonZeroMean(D) == dTotal / (nrow * ncol - 1));
    }

    SECTION("scalar multiplication/division")
    {
        REQUIRE(gaps::algo::sum(gaps::algo::scalarMultiple(v, 3.5))
            == 3.5 * 300.0);

        float vsqSum = 24.0 * 25.0 * (2.0 * 24.0 + 1.0) / 6.0;
        REQUIRE(gaps::algo::sum(gaps::algo::squaredScalarMultiple(v, 4.0))
            == 16.0 * vsqSum);
/*
        REQUIRE(gaps::algo::sum(gaps::algo::scalarDivision(v, 1.3))
            == 300.0 / 1.3);

        REQUIRE(gaps::algo::sum(gaps::algo::squaredScalarDivision(v, 1.3))
            == vsqSum / (1.3 * 1.3));*/
    }
/*        REQUIRE(gaps::algo::sum(D) == 467500);
        
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
*/
}

