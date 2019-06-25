#ifndef __COGAPS_ALPHA_PARAMETERS_H__
#define __COGAPS_ALPHA_PARAMETERS_H__

struct OptionalFloat;
class GapsRng;

struct AlphaParameters
{
    AlphaParameters(float t_s, float t_smu);
    AlphaParameters operator+(const AlphaParameters &other) const;
    AlphaParameters operator*(float v) const;
    void operator*=(float v);

    float s;
    float s_mu;
};

OptionalFloat gibbsMass(AlphaParameters alpha, float a, float b, GapsRng *rng);
OptionalFloat gibbsMass(AlphaParameters alpha, float a, float b, GapsRng *rng, float lambda);

#endif // __COGAPS_ALPHA_PARAMETERS_H__