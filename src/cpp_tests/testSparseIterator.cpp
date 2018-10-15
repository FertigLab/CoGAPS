#include "catch.h"
#include "../math/Random.h"
#include "../math/VectorMath.h"
#include "../data_structures/Matrix.h"
#include "../data_structures/SparseVector.h"
#include "../data_structures/SparseMatrix.h"
#include "../data_structures/HybridMatrix.h"
#include "../data_structures/SparseIterator.h"

TEST_CASE("Test SparseIterator.h - One Dimensional")
{
#ifdef GAPS_INTERNAL_TESTS
    SECTION("Simple Case")
    {
        SparseVector v(10);
        v.insert(0, 1.f);
        v.insert(4, 5.f);
        v.insert(7, 8.f);
        v.insert(9, 10.f);

        SparseIterator it(v);
        REQUIRE(it.getValue() == 1.f);
        it.next();
        REQUIRE(it.getValue() == 5.f);
        it.next();
        REQUIRE(it.getValue() == 8.f);
        it.next();
        REQUIRE(it.getValue() == 10.f);
        it.next();
        REQUIRE(it.atEnd());
    }
#endif

    SECTION("Test Identical Sums")
    {
        GapsRng::setSeed(123);
        GapsRng rng;
        Matrix ref(10, 15);
        for (unsigned i = 0; i < ref.nRow(); ++i)
        {
            for (unsigned j = 0; j < ref.nCol(); ++j)
            {
                ref(i,j) = (i + j) * (rng.uniform() < 0.5f ? 0.f : 1.f);
            }
        }

        SparseMatrix mat(ref, false, false, std::vector<unsigned>());

        for (unsigned j = 0; j < ref.nCol(); ++j)
        {
            float colSum = 0.f;
            SparseIterator it(mat.getCol(j));
            while (!it.atEnd())
            {
                colSum += it.getValue();
                it.next();
            }
            REQUIRE(colSum == gaps::sum(ref.getCol(j)));
        }
    }
}

TEST_CASE("Test SparseIterator.h - Two Dimensional")
{
#ifdef GAPS_INTERNAL_TESTS
   SECTION("Simple Case")
    {
        SparseVector sv(10);
        sv.insert(0, 1.f);
        sv.insert(4, 5.f);
        sv.insert(7, 8.f);
        sv.insert(9, 10.f);

        HybridVector hv(10);
        hv.add(4, 3.f);
        hv.add(5, 4.f);
        hv.add(6, 5.f);
        hv.add(7, 6.f);
        
        SparseIteratorTwo it(sv, hv);
        REQUIRE(it.getValue_1() == 5.f);
        REQUIRE(it.getValue_2() == 3.f);
        it.next();
        REQUIRE(it.getValue_1() == 8.f);
        REQUIRE(it.getValue_2() == 6.f);
        it.next();
        REQUIRE(it.atEnd());
    }
#endif

    // could this fail because of SIMD?
    SECTION("Test Identical Dot Products")
    {
        GapsRng::setSeed(123);
        GapsRng rng;
        Matrix ref(100, 25);
        HybridMatrix hMat(ref.nRow(), ref.nCol());
        for (unsigned i = 0; i < ref.nRow(); ++i)
        {
            for (unsigned j = 0; j < ref.nCol(); ++j)
            {
                ref(i,j) = (i + j) * (rng.uniform() < 0.5f ? 0.f : 1.f);
                hMat.add(i, j, ref(i,j));
            }
        }
        SparseMatrix sMat(ref, false, false, std::vector<unsigned>());

        for (unsigned j1 = 0; j1 < ref.nCol(); ++j1)
        {
            for (unsigned j2 = j1; j2 < ref.nCol(); ++j2)
            {
                float dot = 0.f;
                SparseIteratorTwo it(sMat.getCol(j1), hMat.getCol(j2));
                while (!it.atEnd())
                {
                    dot += it.getValue_1() * it.getValue_2();
                    it.next();
                }
                REQUIRE(dot == gaps::dot(ref.getCol(j1), ref.getCol(j2)));
            }
        }
    }
}

static float tripleProduct(const Vector &v1, const Vector &v2, const Vector &v3)
{
    float prod = 0.f;
    for (unsigned i = 0; i < v1.size(); ++i)
    {
        prod += v1[i] * v2[i] * v3[i];
    }
    return prod;
}   

TEST_CASE("Test SparseIterator.h - Three Dimensional")
{
#ifdef GAPS_INTERNAL_TESTS
   SECTION("Simple Case")
    {
        SparseVector sv(10);
        sv.insert(0, 1.f);
        sv.insert(4, 5.f);
        sv.insert(7, 8.f);
        sv.insert(8, 9.f);
        sv.insert(9, 10.f);

        HybridVector hv1(10);
        hv1.add(4, 3.f);
        hv1.add(5, 4.f);
        hv1.add(6, 5.f);
        hv1.add(7, 6.f);
        hv1.add(9, 7.f);

        HybridVector hv2(10);
        hv2.add(0, 5.f);
        hv2.add(4, 6.f);
        hv2.add(8, 7.f);
        hv2.add(9, 8.f);
        
        SparseIteratorThree it(sv, hv1, hv2);
        REQUIRE(it.getValue_1() == 5.f); // 4
        REQUIRE(it.getValue_2() == 3.f);
        REQUIRE(it.getValue_3() == 6.f);
        it.next();
        REQUIRE(it.getValue_1() == 10.f); // 9
        REQUIRE(it.getValue_2() == 7.f);
        REQUIRE(it.getValue_3() == 8.f);
        it.next();
        REQUIRE(it.atEnd());
    }
#endif

    SECTION("Test Identical Triple Products")
    {
        GapsRng::setSeed(123);
        GapsRng rng;
        Matrix ref(100, 25);
        HybridMatrix hMat(ref.nRow(), ref.nCol());
        for (unsigned i = 0; i < ref.nRow(); ++i)
        {
            for (unsigned j = 0; j < ref.nCol(); ++j)
            {
                ref(i,j) = (i + j) * (rng.uniform() < 0.5f ? 0.f : 1.f);
                hMat.add(i, j, ref(i,j));
            }
        }
        SparseMatrix sMat(ref, false, false, std::vector<unsigned>());

        for (unsigned j1 = 0; j1 < ref.nCol(); ++j1)
        {
            for (unsigned j2 = j1; j2 < ref.nCol(); ++j2)
            {
                for (unsigned j3 = j2; j3 < ref.nCol(); ++j3)
                {
                    float prod = 0.f;
                    SparseIteratorThree it(sMat.getCol(j1), hMat.getCol(j2),
                        hMat.getCol(j3));
                    while (!it.atEnd())
                    {
                        prod += it.getValue_1() * it.getValue_2() * it.getValue_3();
                        it.next();
                    }
                    REQUIRE(prod == tripleProduct(ref.getCol(j1),
                        ref.getCol(j2), ref.getCol(j3)));
                }
            }
        }
    }
}