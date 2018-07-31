#ifndef __COGAPS_RANDOM_H__
#define __COGAPS_RANDOM_H__

#include "../Archive.h"

#include <fstream>
#include <stdint.h>
#include <vector>

namespace gaps
{
    namespace random
    {
        void setSeed(uint32_t seed);

        float uniform();
        float uniform(float a, float b);
        uint64_t uniform64();
        uint64_t uniform64(uint64_t a, uint64_t b);

        float exponential(float lambda);
        float normal(float mean, float var);
        int poisson(float lambda);

        float inverseNormSample(float a, float b, float mean, float sd);
        float inverseGammaSample(float a, float b, float mean, float sd);

        float d_gamma(float d, float shape, float scale);
        float p_gamma(float p, float shape, float scale);
        float q_gamma(float q, float shape, float scale);
        float d_norm(float d, float mean, float sd);
        float q_norm(float q, float mean, float sd);
        float p_norm(float p, float mean, float sd);

        void save(Archive &ar);
        void load(Archive &ar);
    } // namespace random
} // namespace gaps

#endif // __COGAPS_RANDOM_H__
