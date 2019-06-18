#ifndef __COGAPS_RANDOM_H__
#define __COGAPS_RANDOM_H__

#include "../utils/Archive.h"

#include <fstream>
#include <stdint.h>
#include <vector>

class GapsRng;
class Xoroshiro128plus;
class GapsRandomState;

struct OptionalFloat
{
public :

    OptionalFloat();
    OptionalFloat(float f); // NOLINT

    float value();
    bool hasValue() const;

private :

    float mValue;
    bool mHasValue;
};

// TODO allow this to be rolled back
// PCG random number generator
// This is constructed with a seed pulled from the global state
class GapsRng
{
public:

    explicit GapsRng(GapsRandomState *randState);

    float uniform();
    float uniform(float a, float b);

    uint32_t uniform32();
    uint32_t uniform32(uint32_t a, uint32_t b);

    uint64_t uniform64();
    uint64_t uniform64(uint64_t a, uint64_t b);

    int poisson(double lambda);
    float exponential(float lambda);

    OptionalFloat truncNormal(float a, float b, float mean, float sd);
    float truncGammaUpper(float b, float scale); // shape hardcoded to 2

    friend Archive& operator<<(Archive &ar, const GapsRng &gen);
    friend Archive& operator>>(Archive &ar, GapsRng &gen);

#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    const GapsRandomState *mRandState;
    uint64_t mState;

    uint32_t next();
    void advance();
    uint32_t get() const;

    double uniformd();
    int poissonSmall(double lambda);
    int poissonLarge(double lambda);
};

// private class used by GapsRandomState for seeding individual rngs
class Xoroshiro128plus
{
#ifndef GAPS_INTERNAL_TESTS
private:
#else
public:
#endif

    friend class GapsRandomState;

    explicit Xoroshiro128plus(uint64_t seed);

    uint64_t next();
    void rollBackOnce();

    friend Archive& operator<<(Archive &ar, const Xoroshiro128plus &gen);
    friend Archive& operator>>(Archive &ar, Xoroshiro128plus &gen);

    uint64_t mState[2];
    uint64_t mPreviousState[2];
};

// manages random seed and lookup tables for distribution functions, need to
// avoid global variables for multi-threading issues - this random state
// is created at the beginning of execution and passed down to the classes
// that need it
class GapsRandomState
{
public:

    explicit GapsRandomState(unsigned seed);

    uint64_t nextSeed();
    void rollBackOnce();

    // fast distribution calculations using lookup tables
    float p_norm_fast(float p, float mean, float sd) const;
    float q_norm_fast(float q, float mean, float sd) const;

    friend Archive& operator<<(Archive &ar, const GapsRandomState &s);
    friend Archive& operator>>(Archive &ar, GapsRandomState &s);
   
#ifndef GAPS_INTERNAL_TESTS
private:
#endif

    friend class GapsRng;

    Xoroshiro128plus mSeeder;

    float mErfLookupTable[3001];
    float mErfinvLookupTable[5001];
    float mQgammaLookupTable[5001];

    void initLookupTables();

    GapsRandomState(const GapsRandomState&); // = delete (no c++11)
    GapsRandomState& operator=(const GapsRandomState&); // = delete (no c++11)
};

#endif // __COGAPS_RANDOM_H__