#ifndef __COGAPS_RANDOM_H__
#define __COGAPS_RANDOM_H__

#include <stdint.h>
#include <vector>
#include <fstream>

namespace gaps
{

namespace random
{
    void setSeed(uint32_t seed);

    int poisson(double lambda);
    int uniformInt(int a, int b);
    uint64_t uniform64();
    uint64_t uniform64(uint64_t a, uint64_t b);

    double uniform();
    double uniform(double a, double b);
    double normal(double mean, double var);
    double exponential(double lambda);

    double d_gamma(double d, double shape, double scale);
    double p_gamma(double p, double shape, double scale);
    double q_gamma(double q, double shape, double scale);
    double d_norm(double d, double mean, double sd);
    double q_norm(double q, double mean, double sd);
    double p_norm(double p, double mean, double sd);

    void save(std::ofstream &file);
    void load(std::ifstream &file);
}

}

#endif
