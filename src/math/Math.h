#ifndef __COGAPS_MATH_H__
#define __COGAPS_MATH_H__

#include <stdint.h>

namespace gaps
{
    const float epsilon = 1.0e-10f;
    const float pi = 3.14159265358979323846264f;

    float min(float a, float b);
    unsigned min(unsigned a, unsigned b);
    uint64_t min(uint64_t a, uint64_t b);

    float max(float a, float b);
    unsigned max(unsigned a, unsigned b);
    uint64_t max(uint64_t a, uint64_t b);
}

#endif