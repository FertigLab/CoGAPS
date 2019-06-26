#include "Random.h"
#include "Math.h"
#include "../utils/Archive.h"
#include "../utils/GapsAssert.h"

#include <algorithm>
#include <stdint.h>
#include <cmath>
#include <limits>

const float maxU32AsFloat = static_cast<float>(std::numeric_limits<uint32_t>::max());
const double maxU32AsDouble = static_cast<double>(std::numeric_limits<uint32_t>::max());

/////////////////////////////// OptionalFloat //////////////////////////////////

OptionalFloat::OptionalFloat() : mValue(0.f), mHasValue(false) {}

OptionalFloat::OptionalFloat(float f) : mValue(f), mHasValue(true) {}

float OptionalFloat::value()
{
    return mValue;
}

bool OptionalFloat::hasValue() const
{
    return mHasValue;
}

//////////////////////////////// GapsRng ///////////////////////////////////////

GapsRng::GapsRng(GapsRandomState *randState)
:
mRandState(randState),
mState(randState->nextSeed())
{
    advance();
}

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

// inclusive of a and b
uint32_t GapsRng::uniform32(uint32_t a, uint32_t b)
{
    if (b == a)
    {
        return a;
    }
    uint32_t range = b + 1 - a;
    uint32_t x = uniform32();
    uint32_t iPart = std::numeric_limits<uint32_t>::max() / range;
    while (x >= range * iPart)
    {
        x = uniform32();
    }
    uint32_t unif = x / iPart + a;
    GAPS_ASSERT(unif >= a);
    GAPS_ASSERT(unif <= b);
    return unif;
}

uint64_t GapsRng::uniform64()
{
    uint64_t high = (static_cast<uint64_t>(uniform32()) << 32) & 0xFFFFFFFF00000000ull;
    uint64_t low = uniform32();
    return high | low;
}

// inclusive of a and b
uint64_t GapsRng::uniform64(uint64_t a, uint64_t b)
{
    if (b == a)
    {
        return a;
    }
    uint64_t range = b + 1 - a;
    uint64_t x = uniform64();
    uint64_t iPart = std::numeric_limits<uint64_t>::max() / range;
    while (x >= range * iPart)
    {
        x = uniform64();
    }
    uint64_t unif = x / iPart + a;
    GAPS_ASSERT(unif >= a);
    GAPS_ASSERT(unif <= b);
    return unif;
}

int GapsRng::poisson(double lambda)
{
    return lambda <= 5.0 ? poissonSmall(lambda) : poissonLarge(lambda);
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

// fails if too far in tail
OptionalFloat GapsRng::truncNormal(float a, float b, float mean, float sd)
{
    float pLower = mRandState->p_norm_fast(a, mean, sd);
    float pUpper = mRandState->p_norm_fast(b, mean, sd);
    if (!(pLower > 0.95f || pUpper < 0.05f)) // too far in tail
    {
        GAPS_ASSERT(pLower > 0.f);
        GAPS_ASSERT(pUpper < 1.f);
        float z = mRandState->q_norm_fast(uniform(pLower, pUpper), mean, sd); 
        z = gaps::max(a, gaps::min(z, b));
        return z;
    }
    return OptionalFloat();
}

// shape is hardcoded to 2 since it never changes
float GapsRng::truncGammaUpper(float b, float scale)
{
    float upper = 1.f - std::exp(-b / scale) * (1.f + b / scale);
    const unsigned ndx = static_cast<unsigned>(uniform(0.f, upper * 5000.f));
    GAPS_ASSERT(ndx < Q_GAMMA_LOOKUP_TABLE_SIZE);
    return mRandState->mQgammaLookupTable[ndx] * scale;
}

Archive& operator<<(Archive &ar, const GapsRng &gen)
{
    ar << gen.mState;
    return ar;
}

Archive& operator>>(Archive &ar, GapsRng &gen)
{
    ar >> gen.mState;
    return ar;
}

///////////////////////////// Xoroshiro128plus /////////////////////////////////

static uint64_t rotl(const uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

Xoroshiro128plus::Xoroshiro128plus(uint64_t seed)
{
    mState[0] = seed|1;
    mState[1] = seed|1;
    for (unsigned i = 0; i < 5000; ++i) // warmup
    {
        next();
    }
}

uint64_t Xoroshiro128plus::next()
{
    mPreviousState[0] = mState[0];
    mPreviousState[1] = mState[1];
    const uint64_t s0 = mState[0];
    uint64_t s1 = mState[1];
    uint64_t result = s0 + s1;
    s1 ^= s0;
    mState[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
    mState[1] = rotl(s1, 37); // c
    return result;
}

void Xoroshiro128plus::rollBackOnce()
{
    mState[0] = mPreviousState[0];
    mState[1] = mPreviousState[1];
}

Archive& operator<<(Archive &ar, const Xoroshiro128plus &gen)
{
    ar << gen.mState[0] << gen.mState[1];
    return ar;
}

Archive& operator>>(Archive &ar, Xoroshiro128plus &gen)
{
    ar >> gen.mState[0] >> gen.mState[1];
    return ar;
}

///////////////////////////// GapsRandomState //////////////////////////////////

GapsRandomState::GapsRandomState(unsigned seed) : mSeeder(seed)
{
    initLookupTables();
}

void GapsRandomState::initLookupTables()
{
    // erf
    for (unsigned i = 0; i < ERF_LOOKUP_TABLE_SIZE; ++i)
    {
        float x = static_cast<float>(i) / 1000.f;
        mErfLookupTable[i] = 2.f * gaps::p_norm(x * gaps::sqrt2, 0.f, 1.f) - 1.f;
    }
    GAPS_ASSERT(mErfLookupTable[ERF_LOOKUP_TABLE_SIZE - 1] < 1.f);

    // erfinv
    for (unsigned i = 0; i < ERF_INV_LOOKUP_TABLE_SIZE - 1; ++i)
    {
        float x = static_cast<float>(i) / static_cast<float>(ERF_INV_LOOKUP_TABLE_SIZE - 1);
        mErfinvLookupTable[i] = gaps::q_norm((1.f + x) / 2.f, 0.f, 1.f) / gaps::sqrt2;
    }
    mErfinvLookupTable[ERF_INV_LOOKUP_TABLE_SIZE - 1] = gaps::q_norm(1.9998f / 2.f, 0.f, 1.f) / gaps::sqrt2;

    // qgamma
    mQgammaLookupTable[0] = 0.f;
    for (unsigned i = 1; i < Q_GAMMA_LOOKUP_TABLE_SIZE - 1; ++i)
    {
        float x = static_cast<float>(i) / static_cast<float>(Q_GAMMA_LOOKUP_TABLE_SIZE - 1);
        mQgammaLookupTable[i] = gaps::q_gamma(x, 2.f, 1.f);
    }
    mQgammaLookupTable[Q_GAMMA_LOOKUP_TABLE_SIZE - 1] = gaps::q_gamma(0.9998f, 2.f, 1.f);
}

uint64_t GapsRandomState::nextSeed()
{
    return mSeeder.next();
}

void GapsRandomState::rollBackOnce()
{
    mSeeder.rollBackOnce();
}

float GapsRandomState::p_norm_fast(float p, float mean, float sd) const
{
    float term = (p - mean) / (sd * gaps::sqrt2);
    float erf = 0.f;
    if (term < 0.f)
    {
        term = gaps::max(term, -3.f);
        const unsigned ndx = static_cast<unsigned>(-term * 1000.f);
        GAPS_ASSERT(ndx < ERF_LOOKUP_TABLE_SIZE);
        erf = -mErfLookupTable[ndx];
    }
    else
    {
        term = gaps::min(term, 3.f);
        const unsigned ndx = static_cast<unsigned>(term * 1000.f);
        GAPS_ASSERT(ndx < ERF_LOOKUP_TABLE_SIZE);
        erf = mErfLookupTable[ndx];
    }
    return 0.5f * (1.f + erf);
}

float GapsRandomState::q_norm_fast(float q, float mean, float sd) const
{
    float term = 2.f * q - 1.f;
    float erfinv = 0.f;
    if (term < 0.f)
    {
        const unsigned ndx = static_cast<unsigned>(-term * (ERF_INV_LOOKUP_TABLE_SIZE - 1));
        GAPS_ASSERT(ndx < ERF_INV_LOOKUP_TABLE_SIZE);
        erfinv = -mErfinvLookupTable[ndx];
    }
    else
    {
        const unsigned ndx = static_cast<unsigned>(term * (ERF_INV_LOOKUP_TABLE_SIZE - 1));
        GAPS_ASSERT(ndx < ERF_INV_LOOKUP_TABLE_SIZE);
        erfinv = mErfinvLookupTable[ndx];
    }
    return mean + sd * gaps::sqrt2 * erfinv;
}

Archive& operator<<(Archive &ar, const GapsRandomState &s)
{
    ar << s.mSeeder;
    return ar;
}

Archive& operator>>(Archive &ar, GapsRandomState &s)
{
    ar >> s.mSeeder;
    return ar;
}