#include "Algorithms.h"
#include "../data_structures/Matrix.h"
#include "../utils/GapsAssert.h"
#include "SIMD.h"

#include <algorithm>

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

float gaps::algo::deltaLL(unsigned size, const float *D, const float *S,
const float *AP, const float *mat, float delta)
{
    const gaps::simd::packedFloat pDelta(delta), two(2.f);
    gaps::simd::packedFloat d, pMat, pD, pAP, pS, partialSum(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        d = pDelta * pMat;
        partialSum += (d * (two * (pD - pAP) - d)) / (two * pS * pS);
    }
    float fd, delLL = partialSum.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        fd = delta * mat[j];
        delLL += (fd * (2.f * (D[j] - AP[j]) - fd)) / (2.f * S[j] * S[j]);
    }
    return delLL;
}

float gaps::algo::deltaLL(unsigned size, const float *D, const float *S,
const float *AP, const float *mat1, float delta1, const float *mat2,
float delta2)
{
    const gaps::simd::packedFloat pDelta1(delta1), pDelta2(delta2), two(2.f);
    gaps::simd::packedFloat d, pMat1, pMat2, pD, pAP, pS, partialSum(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat1.load(mat1 + i);
        pMat2.load(mat2 + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        d = pDelta1 * pMat1 + pDelta2 * pMat2;
        partialSum += (d * (two * (pD - pAP) - d)) / (two * pS * pS);
    }
    float fd, delLL = partialSum.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        fd = delta1 * mat1[j] + delta2 * mat2[j];
        delLL += (fd * (2.f * (D[j] - AP[j]) - fd)) / (2.f * S[j] * S[j]);
    }
    return delLL;
}

// horribly slow, don't call often
RowMatrix gaps::algo::matrixMultiplication(const ColMatrix &A, const RowMatrix &B)
{
    GAPS_ASSERT_MSG(A.nCol() == B.nRow(), A.nCol() << " " << B.nRow());
    RowMatrix temp(A.nRow(), B.nCol());
    for (unsigned i = 0; i < A.nRow(); ++i)
    {
        for (unsigned j = 0; j < B.nCol(); ++j)
        {
            float sum = 0.0;
            for (unsigned k = 0; k < A.nCol(); ++k)
            {
                sum += A(i,k) * B(k,j);
            }
            temp(i,j) = sum;
        }
    }
    return temp;
}