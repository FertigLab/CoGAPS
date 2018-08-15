// [[Rcpp::depends(BH)]]

#include "Math.h"
#include "Random.h"
#include "../GapsAssert.h"

// TODO remove boost dependency

#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>

#include <boost/random/exponential_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <algorithm>
#include <stdint.h>

#ifdef __GAPS_OPENMP__
    #include <omp.h>
#else
    typedef int omp_int_t;
    inline omp_int_t omp_get_thread_num() { return 0; }
    inline omp_int_t omp_get_max_threads() { return 1; }
#endif

#define Q_GAMMA_THRESHOLD 0.000001f
#define Q_GAMMA_MIN_VALUE 0.0

//typedef boost::random::mt19937 RNGType;
//typedef boost::random::mt11213b RNGType; // should be faster

static Xoroshiro128plus seeder;

static std::vector<GapsRng>& rng()
{
    static std::vector<GapsRng> allRngs(omp_get_max_threads());
    return allRngs;
}

void gaps::random::save(Archive &ar)
{
    for (unsigned i = 0; i < rng().size(); ++i)
    {
        ar << seeder;
    }
}

void gaps::random::load(Archive &ar)
{
    for (unsigned i = 0; i < rng().size(); ++i)
    {
        ar >> seeder;
    }
}

void gaps::random::setSeed(uint32_t seed)
{
    seeder.seed(seed);
}

int gaps::random::poisson(float lambda)
{   
    return rng().at(omp_get_thread_num()).poisson(lambda);
}

float gaps::random::exponential(float lambda)
{
    return rng().at(omp_get_thread_num()).exponential(lambda);
}

float gaps::random::uniform()
{
    return rng().at(omp_get_thread_num()).uniform();
}

float gaps::random::uniform(float a, float b)
{
    return rng().at(omp_get_thread_num()).uniform(a,b);
}

uint64_t gaps::random::uniform64()
{
    return rng().at(omp_get_thread_num()).uniform64();
}

uint64_t gaps::random::uniform64(uint64_t a, uint64_t b)
{
    return rng().at(omp_get_thread_num()).uniform64(a,b);
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
    boost::math::gamma_distribution<> gam(shape, scale);
    return quantile(gam, q);
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

double gaps::lgamma(double x)
{
    return boost::math::lgamma(x);
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

///////////////////////// NEW /////////////////////////////////////////////////

const float maxU32AsFloat = static_cast<float>(std::numeric_limits<uint32_t>::max());
const double maxU32AsDouble = static_cast<double>(std::numeric_limits<uint32_t>::max());

static uint64_t rotl(const uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

void Xoroshiro128plus::seed(uint64_t seed)
{
    mState[0] = seed|1;
    mState[1] = seed|1;
    warmup();
}

uint64_t Xoroshiro128plus::next()
{
    uint64_t result = 0;
    #pragma omp critical(RngCreation)
    {
        const uint64_t s0 = mState[0];
        uint64_t s1 = mState[1];
        result = s0 + s1;
        s1 ^= s0;
        mState[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
        mState[1] = rotl(s1, 37); // c
    }
    return result;
}

void Xoroshiro128plus::warmup()
{
    for (unsigned i = 0; i < 5000; ++i)
    {
        next();
    }
}

Archive& operator<<(Archive &ar, Xoroshiro128plus &gen)
{
    ar << gen.mState[0] << gen.mState[1];
    return ar;
}

Archive& operator>>(Archive &ar, Xoroshiro128plus &gen)
{
    ar >> gen.mState[0] >> gen.mState[1];
    return ar;
}

GapsRng::GapsRng() : mState(seeder.next()) {}

uint32_t GapsRng::next()
{
    advance();
    return get();
}

void GapsRng::advance()
{
    mState = mState * 6364136223846793005ull + (54u|1);
}

uint32_t GapsRng::get() const
{
    uint32_t xorshifted = ((mState >> 18u) ^ mState) >> 27u;
    uint32_t rot = mState >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

double GapsRng::uniformd()
{
    return static_cast<double>(uniform32()) / maxU32AsDouble;
}

float GapsRng::uniform()
{
    return static_cast<float>(uniform32()) / maxU32AsFloat;
}

float GapsRng::uniform(float a, float b)
{
    return uniform() * (b - a) + a;
}

uint32_t GapsRng::uniform32()
{
    return next();
}

uint32_t GapsRng::uniform32(uint32_t a, uint32_t b)
{
    uint32_t range = b - a;
    if (range == 0)
    {
        return a;
    }

    uint32_t x = uniform32();
    uint32_t iPart = std::numeric_limits<uint32_t>::max() / range;
    while (x >= range * iPart)
    {
        x = uniform32();
    }
    return x / iPart + a;
}

uint64_t GapsRng::uniform64()
{
    uint64_t high = (static_cast<uint64_t>(uniform32()) << 32) & 0xFFFFFFFF00000000ull;
    uint64_t low = uniform32();
    return high | low;
}

uint64_t GapsRng::uniform64(uint64_t a, uint64_t b)
{
    uint64_t range = b - a;
    if (range == 0)
    {
        return a;
    }

    uint64_t x = uniform64();
    uint64_t iPart = std::numeric_limits<uint64_t>::max() / range;
    while (x >= range * iPart)
    {
        x = uniform64();
    }
    return x / iPart + a;
}

int GapsRng::poisson(double lambda)
{
    if (lambda <= 5.0)
    {
        return poissonSmall(lambda);
    }
    else
    {
        return poissonLarge(lambda);
    }
}

// lambda <= 5
int GapsRng::poissonSmall(double lambda)
{
    int x = 0;
    double p = uniformd();
    double cutoff = std::exp(-lambda);
    while (p >= cutoff)
    {
        p *= uniformd();
        ++x;
    }
    return x;
}

// lambda > 5
int GapsRng::poissonLarge(double lambda)
{
    double c = 0.767 - 3.36 / lambda;
    double beta = gaps::pi_double / sqrt(3.0 * lambda);
    double alpha = beta * lambda;
    double k = std::log(c) - lambda - std::log(beta);

    while(true)
    {
        double u = uniformd();
        double x = (alpha - std::log((1.0 - u) / u)) / beta;
        double n = floor(x + 0.5);
        if (n < 0.0)
        {
            continue;
        }

        double v = uniformd();
        double y = alpha - beta * x;
        double w = 1.0 + std::exp(y);
        double lhs = y + std::log(v / (w * w));
        double rhs = k + n * std::log(lambda) - gaps::lgamma(n+1);
        if (lhs <= rhs)
        {
            return n;
        }
    }
}

float GapsRng::exponential(float lambda)
{
    return -1.f * std::log(uniform()) / lambda;
}

