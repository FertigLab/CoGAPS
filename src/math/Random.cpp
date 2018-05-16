// [[Rcpp::depends(BH)]]

#include "Random.h"
#include "../GapsAssert.h"

// TODO remove boost dependency
// TODO make concurrent random generation safe
// spawn a random generator for each thread
// - no way to acheive consistent results across different nThreads

#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/gamma.hpp>

#include <stdint.h>
#include <omp.h>

#define Q_GAMMA_THRESHOLD 0.000001f
#define Q_GAMMA_MIN_VALUE 0.0

//typedef boost::random::mt19937 RNGType;
typedef boost::random::mt11213b RNGType; // should be faster

static std::vector<RNGType> rng;

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
    unsigned n = omp_get_max_threads();

    rng.push_back(RNGType());
    rng[0].seed(seed);

    boost::random::uniform_int_distribution<uint32_t> seedDist(0,
        std::numeric_limits<uint32_t>::max());

    for (unsigned i = 1; i < n; ++i)
    {
        rng.push_back(RNGType());
        uint32_t newSeed = seedDist(rng[i-1]);
        rng[i].seed(newSeed);
    }
}

float gaps::random::normal(float mean, float var)
{
    boost::random::normal_distribution<float> dist(mean, var);
    return dist(rng[omp_get_thread_num()]);
}

int gaps::random::poisson(float lambda)
{
    boost::random::poisson_distribution<> dist(lambda);
    return dist(rng[omp_get_thread_num()]);
}

float gaps::random::exponential(float lambda)
{
    boost::random::exponential_distribution<> dist(lambda);
    return dist(rng[omp_get_thread_num()]);
}

float gaps::random::uniform()
{
    float ret = 0.f;
    boost::random::uniform_01<RNGType&> u01_dist(rng[omp_get_thread_num()]);
    return u01_dist();
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
        return dist(rng[omp_get_thread_num()]);
    }
}

uint64_t gaps::random::uniform64()
{
    boost::random::uniform_int_distribution<uint64_t> dist(0,
        std::numeric_limits<uint64_t>::max());
    return dist(rng[omp_get_thread_num()]);
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
        return dist(rng[omp_get_thread_num()]);
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

float gaps::random::inverseNormSample(float a, float b, float mean, float sd)
{
    float u = gaps::random::uniform(a, b);
    while (u == 0.f || u == 1.f)
    {
        u = gaps::random::uniform(a, b);
    }
    return gaps::random::q_norm(u, mean, sd);
}

float gaps::random::inverseGammaSample(float a, float b, float mean, float sd)
{
    float u = gaps::random::uniform(a, b);
    while (u == 0.f || u == 1.f)
    {
        u = gaps::random::uniform(a, b);
    }
    return gaps::random::q_gamma(u, mean, sd);
}