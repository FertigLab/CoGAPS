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

#define Q_GAMMA_THRESHOLD 0.01
#define Q_GAMMA_MIN_VALUE 0.0

//typedef boost::random::mt19937 RNGType;
typedef boost::random::mt11213b RNGType; // should be faster
static RNGType rng;

void Random::setSeed(uint32_t seed)
{
    rng.seed(seed);
}

double Random::normal(double mean, double var)
{
    boost::random::normal_distribution<double> dist(mean, var);
    return dist(rng);
}

double Random::uniform()
{
    boost::random::uniform_01<RNGType&> dist(rng); // could be moved out
    return dist();
}

int Random::poisson(double lambda)
{
    boost::random::poisson_distribution<> dist(lambda);
    return dist(rng);
}

double Random::exponential(double lambda)
{
    boost::random::exponential_distribution<> dist(lambda);
    return dist(rng);
}

uint64_t Random::uniform64()
{
    boost::random::uniform_int_distribution<uint64_t> dist(0,
        std::numeric_limits<uint64_t>::max());
    return dist(rng);
}

/*double Random::p_exp(double p, double rate)
>>>>>>> develop
{
    boost::math::exponential_distribution<> exp(rate);
    return cdf(exp, p);
}

double Random::q_exp(double q, double rate)
{
    boost::math::exponential_distribution<> exp(rate);
    return quantile(exp, q);
}*/

double Random::d_gamma(double d, double shape, double scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return pdf(gam, d);
}

double Random::p_gamma(double p, double shape, double scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return cdf(gam, p);
}

double Random::q_gamma(double q, double shape, double scale)
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

double Random::d_norm(double d, double mean, double sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return pdf(norm, d);
}

double Random::q_norm(double q, double mean, double sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return quantile(norm, q);
}

double Random::p_norm(double p, double mean, double sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return cdf(norm, p);
}

int Random::uniformInt(int a, int b)
{
    if (a > b)
    {
        throw std::invalid_argument("uniformInt: invalid range\n");
    }
    else if (a == b)
    {
        return a;
    }
    else
    {
        boost::random::uniform_int_distribution<> dist(a,b);
        return dist(rng);
    }
}

double Random::uniform(double a, double b)
{
    if (a > b)
    {
        throw std::invalid_argument("uniform: invalid range\n");
    }
    else if (a == b)
    {
        return a;
    }
    else
    {
        boost::random::uniform_real_distribution<> dist(a,b);
        return dist(rng);
    }
}