#include "AlphaParameters.h"

AlphaParameters::AlphaParameters(float t_s, float t_smu)
    : s(t_s), s_mu(t_smu)
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
    s *= v; s_mu *= v;
}