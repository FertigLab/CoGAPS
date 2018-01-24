#include "catch.h"
#include "../Matrix.h"
#include "../Algorithms.h"

#define MAT_SUM(nR, nC) ((nR + nC - 2) * nR * nC / 2.f)

TEST_CASE("Test Algorithms.h")
{
    unsigned nrow = 25;
    unsigned ncol = 20;
    unsigned nfactor = 7;
    Vector v(nrow);
    TwoWayMatrix D(nrow, ncol), S(nrow, ncol), AP(nrow, ncol);
    ColMatrix A(nrow, nfactor);
    RowMatrix P(nfactor, ncol);

    for (unsigned i = 0; i < nrow; ++i)
    {
        v[i] = i;
        for (unsigned j = 0; j < ncol; ++j)
        {
            D.set(i, j, (i + j) / 10.f);
            S.set(i, j, (i + j) / 100.f);
            for (unsigned k = 1; k < nfactor; ++k)
            {
                A(i,k) = i + k;
                P(k,j) = k + j;
            }
        }
    }
    gaps::algo::matrixMultiplication(AP, A, P);

    float Dsum = MAT_SUM(nrow, ncol) / 10.f;
    float Ssum = MAT_SUM(nrow, ncol) / 100.f;

    SECTION("sum")
    {
        REQUIRE(gaps::algo::sum(v) == 300);
        REQUIRE(gaps::algo::sum(D) == Approx(Dsum).epsilon(0.001));
        REQUIRE(gaps::algo::sum(S) == Approx(Ssum).epsilon(0.001));
    }

    SECTION("mean")
    {
        REQUIRE(gaps::algo::mean(D) == Dsum / (nrow * ncol));
        REQUIRE(gaps::algo::nonZeroMean(D) == Dsum / (nrow * ncol - 1));
    }

    SECTION("scalar multiplication/division")
    {
        float vsqSum = 24.f * 25.f * (2.f * 24.f + 1.f) / 6.f;

        REQUIRE(gaps::algo::sum(gaps::algo::scalarMultiple(v, 3.5f))
            == 3.5f * 300.f);
        REQUIRE(gaps::algo::sum(gaps::algo::squaredScalarMultiple(v, 4.f))
            == 16.f * vsqSum);
        REQUIRE(gaps::algo::sum(gaps::algo::scalarDivision(v, 1.3f))
            == 300.f / 1.3f);
        REQUIRE(gaps::algo::sum(gaps::algo::squaredScalarDivision(v, 1.3f))
            == Approx(vsqSum / (1.3f * 1.3f)).epsilon(0.01));
    }

    SECTION("is row/col zero")
    {
        REQUIRE(gaps::algo::isRowZero(P, 0));
        REQUIRE(!gaps::algo::isRowZero(P, 1));
        
        REQUIRE(gaps::algo::isColZero(A, 0));
        REQUIRE(!gaps::algo::isColZero(A, 1));
    }
}