#ifndef __COGAPS_ALPHA_PARAMETERS_H__
#define __COGAPS_ALPHA_PARAMETERS_H__

#include "../math/Random.h"

struct AlphaParameters
{
    float s;
    float s_mu;
    
    AlphaParameters(float t_s, float t_smu);

    AlphaParameters operator+(const AlphaParameters &other) const;
    AlphaParameters operator*(float v) const;
    void operator*=(float v);
};

OptionalFloat gibbsMass(AlphaParameters alpha, float a, float b, GapsRng *rng);
OptionalFloat gibbsMass(AlphaParameters alpha, float a, float b, GapsRng *rng, float lambda);

#endif // __COGAPS_ALPHA_PARAMETERS_H__