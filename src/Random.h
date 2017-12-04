#ifndef __COGAPS_RANDOM_H__
#define __COGAPS_RANDOM_H__

#include <stdint.h>

namespace Random
{
    void setSeed(uint32_t seed);

    int poisson(double lambda);
    int uniformInt(int a, int b);
    uint64_t uniform64();

    double uniform();
    double uniform(double a, double b);
    double normal(double mean, double var);
    double exponential(double lambda);

    double p_exp(double p, double rate);
    double q_exp(double q, double rate);
    double d_gamma(double d, double shape, double scale);
    double p_gamma(double p, double shape, double scale);
    double q_gamma(double q, double shape, double scale);
    double d_norm(double d, double mean, double sd);
    double q_norm(double q, double mean, double sd);
    double p_norm(double p, double mean, double sd);
}

#endif
