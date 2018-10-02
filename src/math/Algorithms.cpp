#include "Algorithms.h"
#include "../data_structures/Matrix.h"
#include "../utils/GapsAssert.h"
#include "SIMD.h"

#include <algorithm>

#define GAPS_SQ(x) ((x) * (x))

////////////////////////////// AlphaParameters /////////////////////////////////

AlphaParameters::AlphaParameters(float inS, float inS_mu)
    : s(inS), s_mu(inS_mu) 
{}

AlphaParameters AlphaParameters::operator+(const AlphaParameters &other) const
{
    return AlphaParameters(s + other.s, s_mu - other.s_mu); // not a typo
}

AlphaParameters AlphaParameters::operator*(float v) const
{
    return AlphaParameters(s * v, s_mu * v);
}

void AlphaParameters::operator*=(float v)
{
    s *= v;
    s_mu *= v;
}

//////////////////////////////// Algorithms ////////////////////////////////////

float gaps::min(float a, float b)
{
    return a < b ? a : b;
}

unsigned gaps::min(unsigned a, unsigned b)
{
    return a < b ? a : b;
}

uint64_t gaps::min(uint64_t a, uint64_t b)
{
    return a < b ? a : b;
}

float gaps::max(float a, float b)
{
    return a < b ? b : a;
}

unsigned gaps::max(unsigned a, unsigned b)
{
    return a < b ? b : a;
}

uint64_t gaps::max(uint64_t a, uint64_t b)
{
    return a < b ? b : a;
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

Vector gaps::algo::elementSq(Vector vec)
{
    unsigned sz = vec.size();

    for (unsigned i = 0; i < sz; ++i)
    {
        vec[i] *= vec[i];
    }
    return vec;
}

float gaps::algo::max(const Vector &vec)
{
    unsigned sz = vec.size();

    float max = vec[0];
    for (unsigned i = 0; i < sz; ++i)
    {
        max = gaps::max(max, vec[i]);
    }
    return max;
}

float gaps::algo::sum(const Vector &vec)
{
    unsigned sz = vec.size();

    float sum = 0.f;
    for (unsigned i = 0; i < sz; ++i)
    {
        sum += vec[i];
    }
    return sum;
}

float gaps::algo::sum(const ColMatrix &mat)
{
    unsigned nc = mat.nCol();

    float sum = 0.f;
    for (unsigned j = 0; j < nc; ++j)
    {
        sum += gaps::algo::sum(mat.getCol(j));
    }
    return sum;
}

float gaps::algo::mean(const ColMatrix &mat)
{
    return gaps::algo::sum(mat) / static_cast<float>(mat.nRow() * mat.nCol());
}

float gaps::algo::nonZeroMean(const ColMatrix &mat)
{
    unsigned nr = mat.nRow();
    unsigned nc = mat.nCol();

    float sum = 0.f;
    unsigned nNonZeros = 0;
    for (unsigned j = 0; j < nc; ++j)
    {
        for (unsigned i = 0; i < nr; ++i)
        {
            if (mat(i,j) != 0.f)
            {
                nNonZeros++;
                sum += mat(i,j);
            }
        }
    }
    return sum / static_cast<float>(nNonZeros);
}

ColMatrix gaps::algo::pmax(ColMatrix mat, float factor)
{
    unsigned nr = mat.nRow();
    unsigned nc = mat.nCol();

    for (unsigned j = 0; j < nc; ++j)
    {
        for (unsigned i = 0; i < nr; ++i)
        {
            mat(i,j) = gaps::max(mat(i,j) * factor, factor);
        }
    }
    return mat;
}

ColMatrix gaps::algo::matrixMultiplication(const ColMatrix &A, const ColMatrix &BT)
{
    GAPS_ASSERT_MSG(A.nCol() == BT.nCol(), A.nCol() << " " << BT.nCol());
    unsigned nr = A.nRow();
    unsigned nc = BT.nRow();
    unsigned nf = A.nCol(); // inner dimension

    ColMatrix temp(nr, nc);
    for (unsigned i = 0; i < nr; ++i)
    {
        for (unsigned j = 0; j < nc; ++j)
        {
            temp(i,j) = 0.f;
            for (unsigned k = 0; k < nf; ++k)
            {
                temp(i,j) += A(i,k) * BT(j,k);
            }
        }
    }
    return temp;
}

void gaps::algo::copyTranspose(ColMatrix *dest, const ColMatrix &src,
unsigned nThreads)
{
    GAPS_ASSERT(dest->nRow() == src.nCol());
    GAPS_ASSERT(dest->nCol() == src.nRow());
    unsigned nc = src.nCol();
    unsigned nr = src.nRow();

    #pragma omp parallel for num_threads(nThreads)
    for (unsigned j = 0; j < nc; ++j)
    {
        for (unsigned i = 0; i < nr; ++i)
        {
            dest->operator()(j,i) = src(i,j);
        }
    }
}

float gaps::algo::loglikelihood(const ColMatrix &D, const ColMatrix &S,
const ColMatrix &AP)
{
    unsigned nr = D.nRow();
    unsigned nc = D.nCol();

    float chi2 = 0.f;
    for (unsigned i = 0; i < nr; ++i)
    {
        for (unsigned j = 0; j < nc; ++j)
        {
            chi2 += GAPS_SQ(D(i,j) - AP(i,j)) / GAPS_SQ(S(i,j));
        }
    }
    return chi2 / 2.f;
}

ColMatrix gaps::algo::computeStdDev(ColMatrix stdMat, const ColMatrix &meanMat,
unsigned nUpdates)
{
    GAPS_ASSERT(nUpdates > 1);
    unsigned nr = stdMat.nRow();
    unsigned nc = stdMat.nCol();

    for (unsigned j = 0; j < nc; ++j)
    {
        for (unsigned i = 0; i < nr; ++i)
        {
            float meanTerm = meanMat(i,j) * meanMat(i,j) / static_cast<float>(nUpdates);
            float numer = gaps::max(0.f, stdMat(i,j) - meanTerm);
            stdMat(i,j) = std::sqrt(numer / (static_cast<float>(nUpdates) - 1.f));
        }
    }
    return stdMat;
}

// vec is a column of either A or P
AlphaParameters gaps::algo::alphaParameters(const Sparsevector &D,
const Sparsevector &vec, const float *A, const float *P, unsigned size)
{
    // initialize
    float s = -1.f * Z_1[column] * beta;
    float su = 0.f;
    for (unsigned i = 0; i < size; ++i)
    {
        su += A[i] * Z_2[column, i];
    }
    su *= -1.f * beta;

    // iterate over common non-zero entries
    Sparsevector it(D, vec);
    while (!it.atEnd())
    {
        float term1 = it.firstValue() / it.secondValue();
        float term2 = term1 * term1 + it.firstValue() * it.firstValue() * beta;
        float term3 = beta * it.firstValue() - alpha * term1 / it.secondValue();
        s += alpha * term2;
        s_mu += alpha * term1 + term3 * gaps::algo::dotProduct(A, P, size);
        it.next();
    }
    return AlphaParameters(s, s_mu);
}

AlphaParameters gaps::algo::alphaParameters(unsigned size, const float *D,
const float *S, const float *AP, const float *mat)
{
    gaps::simd::PackedFloat pMat, pD, pAP, pS;
    gaps::simd::PackedFloat partialS(0.f), partialS_mu(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        gaps::simd::PackedFloat ratio(pMat / (pS * pS));
        partialS += pMat * ratio;
        partialS_mu += ratio * (pD - pAP);
    }

    float s = partialS.scalar(), s_mu = partialS_mu.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        float ratio = mat[j] / (S[j] * S[j]);
        s += mat[j] * ratio;
        s_mu += ratio * (D[j] - AP[j]);
    }
    return AlphaParameters(s, s_mu);
}

AlphaParameters gaps::algo::alphaParameters(unsigned size, const float *D,
const float *S, const float *AP, const float *mat1, const float *mat2)
{
    gaps::simd::PackedFloat pMat1, pMat2, pD, pAP, pS;
    gaps::simd::PackedFloat partialS(0.f), partialS_mu(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat1.load(mat1 + i);
        pMat2.load(mat2 + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        gaps::simd::PackedFloat ratio((pMat1 - pMat2) / (pS * pS));
        partialS += (pMat1 - pMat2) * ratio;
        partialS_mu += ratio * (pD - pAP);
    }

    float s = partialS.scalar(), s_mu = partialS_mu.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        float ratio = (mat1[j] - mat2[j]) / (S[j] * S[j]);
        s += (mat1[j] - mat2[j]) * ratio;
        s_mu += ratio * (D[j] - AP[j]);
    }
    return AlphaParameters(s, s_mu);
}

AlphaParameters gaps::algo::alphaParametersWithChange(unsigned size,
const float *D, const float *S, const float *AP, const float *mat, float ch)
{
    gaps::simd::PackedFloat pCh(ch);
    gaps::simd::PackedFloat pMat, pD, pAP, pS;
    gaps::simd::PackedFloat partialS(0.f), partialS_mu(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        gaps::simd::PackedFloat ratio(pMat / (pS * pS));
        partialS += pMat * ratio;
        partialS_mu += ratio * (pD - (pAP + pCh * pMat));
    }

    float s = partialS.scalar(), s_mu = partialS_mu.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        float ratio = mat[j] / (S[j] * S[j]);
        s += mat[j] * ratio;
        s_mu += ratio * (D[j] - (AP[j] + ch * mat[j]));
    }
    return AlphaParameters(s, s_mu);
}

