// [[Rcpp::depends(BH)]]

#include "Random.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/exponential.hpp>

#include <stdint.h>

#ifdef GAPS_DEBUG
#include <stdexcept>
#endif

#define Q_GAMMA_THRESHOLD 1E-6
#define Q_GAMMA_MIN_VALUE 0.0

typedef boost::random::mt19937 RNGType;
//typedef boost::random::mt11213b RNGType; // should be faster
static RNGType rng;

#ifdef GAPS_DEBUG

static std::vector<char> randGenTypes;
static std::vector<double> randGenValues;

std::vector<char> gaps::random::getTypes()
{
    return randGenTypes;
}

std::vector<double> gaps::random::getValues()
{
    return randGenValues;
}

#endif

void gaps::random::setSeed(uint32_t seed)
{
    rng.seed(seed);
}

double gaps::random::normal(double mean, double var)
{
    boost::random::normal_distribution<double> dist(mean, var);
#ifdef GAPS_DEBUG
    double ret = dist(rng);
    randGenTypes.push_back('N');
    randGenValues.push_back(ret);
    return ret;
#else
    return dist(rng);
#endif
}

int gaps::random::poisson(double lambda)
{
    boost::random::poisson_distribution<> dist(lambda);
#ifdef GAPS_DEBUG
    double ret = dist(rng);
    randGenTypes.push_back('P');
    randGenValues.push_back(ret);
    return ret;
#else
    return dist(rng);
#endif
}

double gaps::random::exponential(double lambda)
{
    boost::random::exponential_distribution<> dist(lambda);
#ifdef GAPS_DEBUG
    double ret = dist(rng);
    randGenTypes.push_back('E');
    randGenValues.push_back(ret);
    return ret;
#else
    return dist(rng);
#endif
}

double gaps::random::uniform()
{
    boost::random::uniform_01<RNGType&> dist(rng); // could be moved out
#ifdef GAPS_DEBUG
    double ret = dist();
    randGenTypes.push_back('U');
    randGenValues.push_back(ret);
    return ret;
#else
    return dist();
#endif
}

double gaps::random::uniform(double a, double b)
{
    if (a == b)
    {
        return a;
    }
    else if (a < b)
    {
        boost::random::uniform_real_distribution<> dist(a,b);
#ifdef GAPS_DEBUG
        double ret = dist(rng);
        randGenTypes.push_back('u');
        randGenValues.push_back(ret);
        return ret;
    }
    else
    {
        throw std::runtime_error("invalid arguments in uniform64()");
    }
#else
        return dist(rng);
    }
#endif
}

uint64_t gaps::random::uniform64()
{
    boost::random::uniform_int_distribution<uint64_t> dist(0,
        std::numeric_limits<uint64_t>::max());
#ifdef GAPS_DEBUG
    double ret = dist(rng);
    randGenTypes.push_back('L');
    randGenValues.push_back(ret);
    return ret;
#else
    return dist(rng);
#endif
}

uint64_t gaps::random::uniform64(uint64_t a, uint64_t b)
{
    if (a == b)
    {
        return a;
    }
    else if (a < b)
    {
        boost::random::uniform_int_distribution<uint64_t> dist(a,b);
#ifdef GAPS_DEBUG
        double ret = dist(rng);
        randGenTypes.push_back('l');
        randGenValues.push_back(ret);
        return ret;
    }
    else
    {
        throw std::runtime_error("invalid arguments in uniform64()");
    }
#else
        return dist(rng);
    }
#endif
}

double gaps::random::d_gamma(double d, double shape, double scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return pdf(gam, d);
}

double gaps::random::p_gamma(double p, double shape, double scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return cdf(gam, p);
}

double gaps::random::q_gamma(double q, double shape, double scale)
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

double gaps::random::d_norm(double d, double mean, double sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return pdf(norm, d);
}

double gaps::random::q_norm(double q, double mean, double sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return quantile(norm, q);
}

double gaps::random::p_norm(double p, double mean, double sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return cdf(norm, p);
}

