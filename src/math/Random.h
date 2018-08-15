#ifndef __COGAPS_RANDOM_H__
#define __COGAPS_RANDOM_H__

#include "../Archive.h"

#include <fstream>
#include <stdint.h>
#include <vector>

struct OptionalFloat
{
    bool hasValue;
    float value;

    OptionalFloat() : hasValue(false), value(0.f) {}
    OptionalFloat(float f) : hasValue(true), value(f) {}
};

namespace gaps
{
    double lgamma(double x);

    float d_gamma(float d, float shape, float scale);
    float p_gamma(float p, float shape, float scale);
    float q_gamma(float q, float shape, float scale);
    float d_norm(float d, float mean, float sd);
    float q_norm(float q, float mean, float sd);
    float p_norm(float p, float mean, float sd);
}

// used for seeding individual rngs
class Xoroshiro128plus
{
public:

    void seed(uint64_t seed);
    uint64_t next();

private:

    uint64_t mState[2];
    void warmup();

    friend Archive& operator<<(Archive &ar, Xoroshiro128plus &gen);
    friend Archive& operator>>(Archive &ar, Xoroshiro128plus &gen);
};

// PCG random number generator
// This is constructed with a seed pulled from the global seeder
class GapsRng
{
public:

    GapsRng();

    float uniform();
    float uniform(float a, float b);

    uint32_t uniform32();
    uint32_t uniform32(uint32_t a, uint32_t b);

    uint64_t uniform64();
    uint64_t uniform64(uint64_t a, uint64_t b);

    int poisson(double lambda);
    float exponential(float lambda);

    float inverseNormSample(float a, float b, float mean, float sd);
    float inverseGammaSample(float a, float b, float mean, float sd);

    static void setSeed(uint64_t seed);
    static void save(Archive &ar);
    static void load(Archive &ar);
  
private:

    uint64_t mState;

    uint32_t next();
    void advance();
    uint32_t get() const;

    double uniformd();
    int poissonSmall(double lambda);
    int poissonLarge(double lambda);

    friend Archive& operator<<(Archive &ar, GapsRng &gen);
    friend Archive& operator>>(Archive &ar, GapsRng &gen);
};

#endif // __COGAPS_RANDOM_H__
