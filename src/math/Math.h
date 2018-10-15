#ifndef __COGAPS_MATH_H__
#define __COGAPS_MATH_H__

#include <stdint.h>
#include <string>
#include <sstream>

namespace gaps
{
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

} // namespace gaps

#endif // __COGAPS_MATH_H__