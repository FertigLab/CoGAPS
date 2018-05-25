#ifndef __COGAPS_ALGORITHMS_H__
#define __COGAPS_ALGORITHMS_H__

#include "../data_structures/Matrix.h"

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
        float rsu = su - other.su; // weird
        return AlphaParameters(rs, rsu);
    }
};

namespace gaps
{
namespace algo
{
    const float epsilon = 1.0e-10f;
    const float pi = 3.14159265358979323846264f;

    float stringToFloat(const std::string &s);

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
    template<class GenericMatrix>
    float sum(const GenericMatrix &mat);

    template<class GenericMatrix>
    float mean(const GenericMatrix &mat);

    template<class GenericMatrix>
    float nonZeroMean(const GenericMatrix &mat);

    template<class GenericMatrix>
    GenericMatrix computeStdDev(const GenericMatrix &stdMat,
        const GenericMatrix &meanMat, unsigned nUpdates);

    // specific matrix algorithms
    RowMatrix matrixMultiplication(const ColMatrix &A, const RowMatrix &B);

    // chiSq / 2
    template <class Matrix>
    float loglikelihood(const Matrix &D, const Matrix &S,
        const Matrix &AP);

    AlphaParameters alphaParameters(unsigned size, const float *D,
        const float *S, const float *AP, const float *mat);

    AlphaParameters alphaParameters(unsigned size, const float *D,
        const float *S, const float *AP, const float *mat1, const float *mat2);

    float deltaLL(unsigned size, const float *D, const float *S,
        const float *AP, const float *mat, float delta);

    float deltaLL(unsigned size, const float *D, const float *S,
        const float *AP, const float *mat1, float delta1, const float *mat2,
        float delta2);

} // namespace algo
} // namespace gaps

template<class GenericMatrix>
float gaps::algo::sum(const GenericMatrix &mat)
{
    float sum = 0.f;
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            sum += mat(i,j);
        }
    }
    return sum;
}

template<class GenericMatrix>
float gaps::algo::mean(const GenericMatrix &mat)
{
    return gaps::algo::sum(mat) / (mat.nRow() * mat.nCol());
}

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
    return sum / (float)nNonZeros;
}

template<class GenericMatrix>
GenericMatrix gaps::algo::computeStdDev(const GenericMatrix &stdMat,
const GenericMatrix &meanMat, unsigned nUpdates)
{
    GenericMatrix retMat(stdMat.nRow(), stdMat.nCol());
    float meanTerm = 0.f;
    for (unsigned r = 0; r < retMat.nRow(); ++r)
    {
        for (unsigned c = 0; c < retMat.nCol(); ++c)
        {
            meanTerm = meanMat(r,c) * meanMat(r,c) / (float)nUpdates;
            retMat(r,c) = std::sqrt((stdMat(r,c) - meanTerm) / ((float)nUpdates - 1.f));
        }
    }
    return retMat;
}

template <class Matrix>
float gaps::algo::loglikelihood(const Matrix &D, const Matrix &S,
const Matrix &AP)
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