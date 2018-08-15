#ifndef __COGAPS_MATH_H__
#define __COGAPS_MATH_H__

#include <stdint.h>
#include <string>
#include <sstream>

namespace gaps
{
    const float epsilon = 1.0e-10f;
    const float pi = 3.14159265358979323846264f;
    const float pi_double = 3.14159265358979323846264;

    float min(float a, float b);
    unsigned min(unsigned a, unsigned b);
    uint64_t min(uint64_t a, uint64_t b);

    float max(float a, float b);
    unsigned max(unsigned a, unsigned b);
    uint64_t max(uint64_t a, uint64_t b);

    template <class T>
    std::string to_string(T a);
}

template <class T>
std::string gaps::to_string(T a)
{
    std::stringstream ss;
    ss << a;
    return ss.str();
}

#endif