// [[Rcpp::depends(BH)]]

#include "Random.h"
#include <boost/random.hpp>
#include <stdint.h>

//typedef boost::random::mt19937 RNGType;
typedef boost::random::mt11213b RNGType; // should be faster
static RNGType rng;

void Random::setSeed(uint32_t newSeed)
{
    rng.seed(newSeed);
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

double Random::poisson(double lambda)
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