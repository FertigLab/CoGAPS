#include "Algorithms.h"
#include "../data_structures/Matrix.h"
#include "../utils/GapsAssert.h"
#include "SIMD.h"

#include <algorithm>

ColMatrix gaps::algo::pmax(const ColMatrix &mat, float factor)
{
    ColMatrix temp(mat.nRow(), mat.nCol());
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            temp(i,j) = gaps::max(mat(i,j) * factor, factor);
        }
    }
    return temp;
}

float gaps::algo::sum(const ColMatrix &mat, bool transposeOrder)
{
    float sum = 0.f;
    unsigned outer = transposeOrder ? mat.nCol() : mat.nRow();
    unsigned inner = transposeOrder ? mat.nRow() : mat.nCol();
    for (unsigned i = 0; i < outer; ++i)
    {
        for (unsigned j = 0; j < inner; ++j)
        {
            sum += transposeOrder ? mat(j,i) : mat(i,j);
        }
    }
    return sum;
}

float gaps::algo::mean(const ColMatrix &mat, bool transposeOrder)
{
    return gaps::algo::sum(mat, transposeOrder) / (mat.nRow() * mat.nCol());
}

void gaps::algo::copyTranspose(ColMatrix *dest, const ColMatrix &src)
{
    GAPS_ASSERT(dest->nRow() == src.nCol());
    GAPS_ASSERT(dest->nCol() == src.nRow());

    for (unsigned j = 0; j < src.nCol(); ++j)
    {
        for (unsigned i = 0; i < src.nRow(); ++i)
        {
            dest->operator()(j,i) = src(i,j); // TODO test which order is better
        }
    }
}

float gaps::algo::sum(const Vector &vec)
{
    float sum = 0.f;
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        sum += vec[i];
    }
    return sum;
}

float gaps::algo::min(const Vector &vec)
{
    float min = vec[0];
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        min = gaps::min(min, vec[i]);
    }
    return min;
}

float gaps::algo::max(const Vector &vec)
{
    float max = vec[0];
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        max = gaps::max(max, vec[i]);
    }
    return max;
}

float gaps::algo::dot(const Vector &A, const Vector &B)
{
    float dotProd = 0.f;
    for (unsigned i = 0; i < A.size(); ++i)
    {
        dotProd += A[i] * B[i];
    }
    return dotProd;
}

unsigned gaps::algo::whichMin(const Vector &vec)
{
    float min = vec[0];
    unsigned minNdx = 0;
    for (unsigned i = 1; i < vec.size(); ++i)
    {
        if (vec[i] < min)
        {
            min = vec[i];
            minNdx = i;
        }
    }
    return minNdx;
}

Vector gaps::algo::rank(Vector vec)
{
    std::vector< std::pair<float, float> > sortVec(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        sortVec[i] = std::pair<float, float>(vec[i], i);
    }
    
    std::sort(sortVec.begin(), sortVec.end());
    Vector ranks(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        ranks[i] = sortVec[i].second;
    }
    return ranks;
}

Vector gaps::algo::elementSq(Vector vec)
{
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        vec[i] *= vec[i];
    }
    return vec;
}

bool gaps::algo::isVectorZero(const float *vec, unsigned size)
{
    for (unsigned i = 0; i < size; ++i)
    {
        if (vec[i] != 0.f)
        {
            return false;
        }
    }
    return true;
}
    
AlphaParameters gaps::algo::alphaParameters(unsigned size, const float *D,
const float *S, const float *AP, const float *mat)
{
    gaps::simd::packedFloat ratio, pMat, pD, pAP, pS;
    gaps::simd::packedFloat partialS(0.f), partialSU(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        ratio = pMat / pS;
        partialS += ratio * ratio;
        partialSU += (ratio * (pD - pAP)) / pS;
    }
    float fratio, s = partialS.scalar(), su = partialSU.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        fratio = mat[j] / S[j]; // can save one division here by dividing by S^2
        s += fratio * fratio;
        su += (fratio * (D[j] - AP[j])) / S[j];
    }
    return AlphaParameters(s,su);
}

//
AlphaParameters gaps::algo::alphaParameters(unsigned size, const float *D,
const float *S, const float *AP, const float *mat1, const float *mat2)
{
    gaps::simd::packedFloat ratio, pMat1, pMat2, pD, pAP, pS;
    gaps::simd::packedFloat partialS(0.f), partialSU(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat1.load(mat1 + i);
        pMat2.load(mat2 + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        ratio = (pMat1 - pMat2) / pS;
        partialS += ratio * ratio;
        partialSU += ratio * (pD - pAP) / pS;
    }

    float fratio, s = partialS.scalar(), su = partialSU.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        fratio = (mat1[j] - mat2[j]) / S[j];
        s += fratio * fratio;
        su += fratio * (D[j] - AP[j]) / S[j];
    }
    return AlphaParameters(s,su);
}

// horribly slow, don't call often
ColMatrix gaps::algo::matrixMultiplication(const ColMatrix &A, const ColMatrix &BT)
{
    GAPS_ASSERT_MSG(A.nCol() == BT.nCol(), A.nCol() << " " << BT.nCol());
    ColMatrix temp(A.nRow(), BT.nRow());
    for (unsigned i = 0; i < A.nRow(); ++i)
    {
        for (unsigned j = 0; j < BT.nRow(); ++j)
        {
            float sum = 0.0;
            for (unsigned k = 0; k < A.nCol(); ++k)
            {
                sum += A(i,k) * BT(j,k);
            }
            temp(i,j) = sum;
        }
    }
    return temp;
}