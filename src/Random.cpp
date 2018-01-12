// [[Rcpp::depends(BH)]]

#include "Random.h"

#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>

#include <boost/random/mersenne_twister.hpp>

// need -O0 to run in valgrind
#pragma GCC push_options
#pragma GCC optimize ("O0")
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/gamma.hpp>
#pragma GCC pop_options

#include <stdint.h>

#define Q_GAMMA_THRESHOLD 1E-6
#define Q_GAMMA_MIN_VALUE 0.0

typedef boost::random::mt19937 RNGType;
//typedef boost::random::mt11213b RNGType; // should be faster

static RNGType rng;

void gaps::random::save(Archive &ar)
{
    ar << rng;
}

void gaps::random::load(Archive &ar)
{
    ar >> rng;
}

void gaps::random::setSeed(uint32_t seed)
{
    rng.seed(seed);
}

float gaps::random::normal(float mean, float var)
{
    boost::random::normal_distribution<float> dist(mean, var);
    return dist(rng);
}

int gaps::random::poisson(float lambda)
{
    boost::random::poisson_distribution<> dist(lambda);
    return dist(rng);
}

float gaps::random::exponential(float lambda)
{
    boost::random::exponential_distribution<> dist(lambda);
    return dist(rng);
}

float gaps::random::uniform()
{
    boost::random::uniform_01<RNGType&> dist(rng); // could be moved out
    return dist();
}

float gaps::random::uniform(float a, float b)
{
    if (a == b)
    {
        return a;
    }
    else
    {
        boost::random::uniform_real_distribution<> dist(a,b);
        return dist(rng);
    }
}

uint64_t gaps::random::uniform64()
{
    boost::random::uniform_int_distribution<uint64_t> dist(0,
        std::numeric_limits<uint64_t>::max());
    return dist(rng);
}

uint64_t gaps::random::uniform64(uint64_t a, uint64_t b)
{
    if (a == b)
    {
        return a;
    }
    else
    {
        boost::random::uniform_int_distribution<uint64_t> dist(a,b);
        return dist(rng);
    }
}

float gaps::random::d_gamma(float d, float shape, float scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return pdf(gam, d);
}

float gaps::random::p_gamma(float p, float shape, float scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return cdf(gam, p);
}

float gaps::random::q_gamma(float q, float shape, float scale)
{
    if (q < Q_GAMMA_THRESHOLD)
    {
        return Q_GAMMA_MIN_VALUE;
    }
    else
    {
        boost::math::gamma_distribution<> gam(shape, scale);
        return quantile(gam, q);
    }
}

float gaps::random::d_norm(float d, float mean, float sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return pdf(norm, d);
}

float gaps::random::q_norm(float q, float mean, float sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return quantile(norm, q);
}

float gaps::random::p_norm(float p, float mean, float sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return cdf(norm, p);
}
