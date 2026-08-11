#include <testthat.h>
#include "../testthat-tweak.h"
#include "../data_structures/Matrix.h"
#include "../data_structures/SparseMatrix.h"
#include "../data_structures/Vector.h"
#include "../math/MatrixMath.h"
#include "../math/VectorMath.h"

// Helpers that the sampler and GapsStatistics rely on but that had no direct
// coverage: elementSq, dot_diff, mean and sparsity. dot_diff in particular runs
// a SIMD loop, so it is exercised on sizes that are not a multiple of SIMD_INC.

static std::vector<unsigned> sequentialVector(unsigned n)
{
    std::vector<unsigned> vec;
    for (unsigned i = 1; i <= n; ++i) // mimic R indices
    {
        vec.push_back(i);
    }
    return vec;
}

TEST_CASE("gaps::elementSq squares each element","[vectormath][elementsq]")
{
    SECTION("known values")
    {
        Vector v(4);
        v[0] = 0.f; v[1] = 2.f; v[2] = -3.f; v[3] = 1.5f;
        Vector sq(gaps::elementSq(v));

        REQUIRE(sq.size() == v.size());
        REQUIRE(sq[0] == 0.f);
        REQUIRE(sq[1] == 4.f);
        REQUIRE(sq[2] == 9.f);
        REQUIRE(sq[3] == Approx(2.25f));
        // the source vector is not modified
        REQUIRE(v[2] == -3.f);
    }

    SECTION("empty vector")
    {
        Vector v(0);
        Vector sq(gaps::elementSq(v));
        REQUIRE(sq.size() == 0);
    }
}

TEST_CASE("gaps::dot_diff computes sum(a * (b - c))","[vectormath][dotdiff]")
{
    SECTION("known values")
    {
        Vector a(3), b(3), c(3);
        a[0] = 1.f; a[1] = 2.f; a[2] = 3.f;
        b[0] = 4.f; b[1] = 5.f; b[2] = 6.f;
        c[0] = 1.f; c[1] = 1.f; c[2] = 1.f;
        // 1*(4-1) + 2*(5-1) + 3*(6-1) = 3 + 8 + 15 = 26
        REQUIRE(gaps::dot_diff(a, b, c) == Approx(26.f));
    }

    SECTION("agrees with dot(a, b-c) on a size that is not a SIMD multiple")
    {
        const unsigned n = 13;
        Vector a(n), b(n), c(n), diff(n);
        for (unsigned i = 0; i < n; ++i)
        {
            a[i] = static_cast<float>(i) * 0.5f + 1.f;
            b[i] = static_cast<float>(n - i);
            c[i] = static_cast<float>(i) * 0.25f;
            diff[i] = b[i] - c[i];
        }
        REQUIRE(gaps::dot_diff(a, b, c) == Approx(gaps::dot(a, diff)));
    }

    SECTION("b == c gives zero")
    {
        Vector a(5), b(5);
        for (unsigned i = 0; i < 5; ++i)
        {
            a[i] = static_cast<float>(i + 1);
            b[i] = static_cast<float>(i * 3);
        }
        REQUIRE(gaps::dot_diff(a, b, b) == Approx(0.f));
    }
}

TEST_CASE("gaps::mean and gaps::sparsity on a Matrix","[matrixmath][meansparsity]")
{
    // 2x2 with one zero: sum = 6, mean = 6/4, sparsity = 1 - 3/4
    Matrix mat(2, 2);
    mat(0,0) = 1.f; mat(0,1) = 2.f;
    mat(1,0) = 3.f; mat(1,1) = 0.f;

    SECTION("mean is the sum over the number of cells")
    {
        REQUIRE(gaps::sum(mat) == Approx(6.f));
        REQUIRE(gaps::mean(mat) == Approx(1.5f));
    }

    SECTION("sparsity is the fraction of zero cells")
    {
        REQUIRE(gaps::sparsity(mat) == Approx(0.25f));
    }

    SECTION("an all-zero matrix is fully sparse and has mean 0")
    {
        Matrix zero(3, 3);
        REQUIRE(gaps::mean(zero) == Approx(0.f));
        REQUIRE(gaps::sparsity(zero) == Approx(1.f));
    }

    SECTION("a matrix with no zeros has sparsity 0")
    {
        Matrix full(2, 2);
        full(0,0) = 1.f; full(0,1) = 1.f;
        full(1,0) = 1.f; full(1,1) = 1.f;
        REQUIRE(gaps::sparsity(full) == Approx(0.f));
        REQUIRE(gaps::mean(full) == Approx(1.f));
    }
}

TEST_CASE("gaps::sparsity agrees between Matrix and SparseMatrix","[matrixmath][sparsityagree]")
{
    Matrix mat(4, 4);
    for (unsigned j = 0; j < 4; ++j)
    {
        for (unsigned i = 0; i < 4; ++i)
        {
            // leave the diagonal zero -> 4 zeros out of 16
            mat(i,j) = (i == j) ? 0.f : static_cast<float>(i + j + 1);
        }
    }
    SparseMatrix smat(mat, false, false, sequentialVector(0));

    REQUIRE(gaps::sparsity(mat) == Approx(0.25f));
    REQUIRE(gaps::sparsity(smat) == Approx(gaps::sparsity(mat)));
}
