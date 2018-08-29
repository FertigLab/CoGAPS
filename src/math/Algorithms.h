#ifndef __COGAPS_ALGORITHMS_H__
#define __COGAPS_ALGORITHMS_H__

#include "../data_structures/Matrix.h"
#include "Math.h"

#include <cmath>

#define GAPS_SQ(x) ((x) * (x))

struct AlphaParameters
{
    float s;
    float su;
    
    AlphaParameters(float inS, float inSU)
        : s(inS), su(inSU)
    {}

    AlphaParameters operator+(const AlphaParameters &other) const
    {
        float rs = s + other.s;
        float rsu = su - other.su; // not a typo
        return AlphaParameters(rs, rsu);
    }

    void operator*=(float v)
    {
        s *= v;
        su *= v;
    }
};

namespace gaps
{
namespace algo
{
    bool isVectorZero(const float *vec, unsigned size);
    
    // vector algorithms    
    unsigned whichMin(const Vector &vec);
    float sum(const Vector &vec);
    float min(const Vector &vec);
    float max(const Vector &vec);
    float dot(const Vector &A, const Vector &B);
    Vector rank(Vector vec);
    Vector elementSq(Vector vec);
    
    // generic matrix algorithms
    float sum(const ColMatrix &mat, bool transposeOrder=false);
    float mean(const ColMatrix &mat, bool transposeOrder=false);

    template<class GenericMatrix>
    float nonZeroMean(const GenericMatrix &mat);

    template<class GenericMatrix>
    GenericMatrix computeStdDev(const GenericMatrix &stdMat,
        const GenericMatrix &meanMat, unsigned nUpdates);

    ColMatrix pmax(const ColMatrix &mat, float factor);

    // specific matrix algorithms
    ColMatrix matrixMultiplication(const ColMatrix &A, const ColMatrix &B);

    void copyTranspose(ColMatrix *dest, const ColMatrix &src);

    // chiSq / 2
    template <class MatA, class MatB, class MatC>
    float loglikelihood(const MatA &D, const MatB &S,
        const MatC &AP);

    AlphaParameters alphaParameters(unsigned size, const float *D,
        const float *S, const float *AP, const float *mat);

    AlphaParameters alphaParameters(unsigned size, const float *D,
        const float *S, const float *AP, const float *mat1, const float *mat2);

} // namespace algo
} // namespace gaps

template<class GenericMatrix>
float gaps::algo::nonZeroMean(const GenericMatrix &mat)
{
    float sum = 0.f;
    unsigned nNonZeros = 0;
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
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

template<class GenericMatrix>
GenericMatrix gaps::algo::computeStdDev(const GenericMatrix &stdMat,
const GenericMatrix &meanMat, unsigned nUpdates)
{
    GAPS_ASSERT(nUpdates > 1);
    GenericMatrix retMat(stdMat.nRow(), stdMat.nCol());
    for (unsigned r = 0; r < retMat.nRow(); ++r)
    {
        for (unsigned c = 0; c < retMat.nCol(); ++c)
        {
            float meanTerm = meanMat(r,c) * meanMat(r,c) / static_cast<float>(nUpdates);
            float numer = gaps::max(0.f, stdMat(r,c) - meanTerm);
            retMat(r,c) = std::sqrt(numer / (static_cast<float>(nUpdates) - 1.f));
        }
    }
    return retMat;
}

template <class MatA, class MatB, class MatC>
float gaps::algo::loglikelihood(const MatA &D, const MatB &S,
const MatC &AP)
{
    float chi2 = 0.f;
    for (unsigned i = 0; i < D.nRow(); ++i)
    {
        for (unsigned j = 0; j < D.nCol(); ++j)
        {
            chi2 += GAPS_SQ(D(i,j) - AP(i,j)) / GAPS_SQ(S(i,j));
        }
    }
    return chi2 / 2.f;
}

#endif