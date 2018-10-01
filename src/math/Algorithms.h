#ifndef __COGAPS_ALGORITHMS_H__
#define __COGAPS_ALGORITHMS_H__

#include "../data_structures/Matrix.h"

#include <string>
#include <sstream>

struct AlphaParameters
{
    float s;
    float s_mu;
    
    AlphaParameters(float inS, float inS_mu);

    AlphaParameters operator+(const AlphaParameters &other) const;
    AlphaParameters operator*(float v) const;
    void operator*=(float v);
};

namespace gaps
{
    template <class T>
    std::string to_string(T a);

    const float epsilon = 1.0e-5f;
    const float pi = 3.1415926535897932384626433832795f;
    const float pi_double = 3.1415926535897932384626433832795;
    const float sqrt2 = 1.4142135623730950488016887242097f;

    float min(float a, float b);
    unsigned min(unsigned a, unsigned b);
    uint64_t min(uint64_t a, uint64_t b);

    float max(float a, float b);
    unsigned max(unsigned a, unsigned b);
    uint64_t max(uint64_t a, uint64_t b);

namespace algo
{
    // vector algorithms
    bool isVectorZero(const float *vec, unsigned size);
    Vector elementSq(Vector vec);
    float max(const Vector &vec);
    float sum(const Vector &vec);

    // matrix algorithms
    float sum(const ColMatrix &mat);
    float mean(const ColMatrix &mat);
    float nonZeroMean(const ColMatrix &mat);
    ColMatrix pmax(ColMatrix mat, float factor);
    ColMatrix matrixMultiplication(const ColMatrix &A, const ColMatrix &BT);
    void copyTranspose(ColMatrix *dest, const ColMatrix &src, unsigned nThreads=1);
    ColMatrix computeStdDev(ColMatrix stdMat, const ColMatrix &meanMat,
        unsigned nUpdates);

    // Cogaps Algorithms
    float loglikelihood(const ColMatrix &D, const ColMatrix &S,
        const ColMatrix &AP);

    AlphaParameters alphaParameters(unsigned size, const float *D,
        const float *S, const float *AP, const float *mat);

    AlphaParameters alphaParameters(unsigned size, const float *D,
        const float *S, const float *AP, const float *mat1, const float *mat2);

    AlphaParameters alphaParametersWithChange(unsigned size, const float *D,
        const float *S, const float *AP, const float *mat, float ch);

} // namespace algo
} // namespace gaps

template <class T>
std::string gaps::to_string(T a)
{
    std::stringstream ss;
    ss << a;
    return ss.str();
}

#endif // __COGAPS_ALGORITHMS_H__