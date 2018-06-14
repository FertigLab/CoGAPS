#include "catch.h"
#include "../data_structures/Matrix.h"
#include "../math/Algorithms.h"

#define MAT_SUM(nR, nC) ((nR + nC - 2) * nR * nC / 2.f)

TEST_CASE("Test Algorithms.h")
{
    unsigned nrow = 25;
    unsigned ncol = 20;
    unsigned nfactor = 7;
    ColMatrix A(nrow, nfactor);
    RowMatrix P(nfactor, ncol);

    for (unsigned i = 0; i < nrow; ++i)
    {
        for (unsigned j = 0; j < ncol; ++j)
        {
            for (unsigned k = 1; k < nfactor; ++k)
            {
                A(i,k) = i + k;
                P(k,j) = k + j;
            }
        }
    }

    SECTION("is row/col zero")
    {
        REQUIRE(gaps::algo::isVectorZero(P.rowPtr(0), P.nCol()));
        REQUIRE(!gaps::algo::isVectorZero(P.rowPtr(1), P.nCol()));
        
        REQUIRE(gaps::algo::isVectorZero(A.colPtr(0), A.nRow()));
        REQUIRE(!gaps::algo::isVectorZero(A.colPtr(1), A.nRow()));
    }
}
