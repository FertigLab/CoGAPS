#include "catch.h"
#include "../data_structures/HybridMatrix.h"
#include "../data_structures/Matrix.h"
#include "../math/VectorMath.h"
#include "../math/MatrixMath.h"

TEST_CASE("Test HybridMatrix.h")
{
    HybridMatrix mat(100, 250);
    REQUIRE(mat.nRow() == 100);
    REQUIRE(mat.nCol() == 250);

    REQUIRE(gaps::sum(mat.getRow(53)) == 0.f);
    mat.add(53, 100, 74.5f);
    REQUIRE(gaps::sum(mat.getRow(53)) == 74.5f);

    Matrix ref(mat.nRow(), mat.nCol());
    for (unsigned i = 0; i < ref.nRow(); ++i)
    {
        for (unsigned j = 0; j < ref.nCol(); ++j)
        {
            ref(i,j) = i + j;
        }
    }
    mat = ref;

    for (unsigned j = 0; j < mat.nCol(); ++j)
    {
        for (unsigned i = 0; i < mat.nRow(); ++i)
        {
            REQUIRE(mat(i,j) == ref(i,j));
        }
    }
}