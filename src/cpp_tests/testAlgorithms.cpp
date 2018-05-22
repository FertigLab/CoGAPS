#include "catch.h"
#include "../data_structures/Matrix.h"
#include "../math/Algorithms.h"

#define MAT_SUM(nR, nC) ((nR + nC - 2) * nR * nC / 2.f)

#if 0

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

    SECTION("is row/col zero")
    {
        REQUIRE(gaps::algo::isRowZero(P, 0));
        REQUIRE(!gaps::algo::isRowZero(P, 1));
        
        REQUIRE(gaps::algo::isColZero(A, 0));
        REQUIRE(!gaps::algo::isColZero(A, 1));
    }
}

#endif