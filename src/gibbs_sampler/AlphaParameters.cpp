#include "AlphaParameters.h"
#include "../math/Math.h"
#include "../math/Random.h"

#include <cmath>

AlphaParameters::AlphaParameters(float t_s, float t_smu)
    : s(t_s), s_mu(t_smu)
{}

AlphaParameters AlphaParameters::operator+(const AlphaParameters &other) const
{
    return AlphaParameters(s + other.s, s_mu - other.s_mu); // minus sign not a typo
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

OptionalFloat gibbsMass(AlphaParameters alpha, float a, float b, GapsRng *rng)
{
    if (alpha.s > gaps::epsilon)
    {
        float mean = alpha.s_mu / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        return rng->truncNormal(a, b, mean, sd);
    }
    return OptionalFloat();
}

OptionalFloat gibbsMass(AlphaParameters alpha, float a, float b, GapsRng *rng,
float lambda)
{
    if (alpha.s > gaps::epsilon)
    {
        float mean = (alpha.s_mu - lambda) / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        return rng->truncNormal(a, b, mean, sd);
    }
    return OptionalFloat();
}