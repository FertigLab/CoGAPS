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

    double pexp(double p, double rate, bool boolpara1, bool boolpara2);
    double qexp(double q, double rate, bool boolpara1, bool boolpara2);
    double dgamma(double newMass, double shape, double scale, bool boolpara);
    double pgamma(double p, double shape, double scale, bool boolpara1, bool boolpara2);
    double qgamma(double q, double shape, double scale, bool boolpara1, bool boolpara2);
    double dnorm(double u, double mean, double sd, bool unknown);
    double qnorm(double u, double mean, double sd, double INF_Ref, double unknown);
    double pnorm(double u, double mean, double sd, double INF_Ref, double unknown);
}

#endif
