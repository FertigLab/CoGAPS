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

    //double pexp(double p, double rate);
    //double qexp(double q, double rate);
    double dgamma(double d, double shape, double scale);
    double pgamma(double p, double shape, double scale);
    double qgamma(double q, double shape, double scale);
    double dnorm(double d, double mean, double sd);
    double qnorm(double q, double mean, double sd);
    double pnorm(double p, double mean, double sd);
}

#endif
