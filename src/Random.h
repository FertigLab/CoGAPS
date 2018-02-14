#ifndef __COGAPS_RANDOM_H__
#define __COGAPS_RANDOM_H__

#include "Archive.h"

#include <stdint.h>
#include <vector>
#include <fstream>

namespace gaps
{

// rename to math for dist functions
namespace random
{
    void setSeed(uint32_t seed);

    int poisson(float lambda);
    uint64_t uniform64();
    uint64_t uniform64(uint64_t a, uint64_t b);

    float uniform();
    float uniform(float a, float b);
    float normal(float mean, float var);
    float exponential(float lambda);

    float d_gamma(float d, float shape, float scale);
    float p_gamma(float p, float shape, float scale);
    float q_gamma(float q, float shape, float scale);
    float d_norm(float d, float mean, float sd);
    float q_norm(float q, float mean, float sd);
    float p_norm(float p, float mean, float sd);

    void save(Archive &ar);
    void load(Archive &ar);
}

}

#endif
