#ifndef _RANDGEN_H_
#define _RANDGEN_H_

// wrapper class for random number generator

#include <stdint.h>

namespace Random
{
    void setSeed(uint32_t);

    int uniformInt(int, int);
    uint64_t uniform64();

    double uniform();
    double uniform(double, double);
    double normal(double, double);
    double poisson(double);
    double exponential(double);
}

#endif
