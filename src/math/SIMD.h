#ifndef __COGAPS_SIMD_H__
#define __COGAPS_SIMD_H__

#ifndef SSE_INSTR_SET
    #if defined ( __AVX2__ )
        #define SSE_INSTR_SET 8
    #elif defined ( __AVX__ )
        #define SSE_INSTR_SET 7
    #elif defined ( __SSE4_2__ )
        #define SSE_INSTR_SET 6
    #elif defined ( __SSE4_1__ )
        #define SSE_INSTR_SET 5
    #else
        #define SSE_INSTR_SET 0
    #endif
#endif

#if SSE_INSTR_SET > 6
    #define __GAPS_AVX__
    #include <immintrin.h>
#elif SSE_INSTR_SET == 6 || SSE_INSTR_SET == 5
    #define __GAPS_SSE__
    #include <nmmintrin.h>
#endif

namespace gaps
{
namespace simd
{

#if defined( __GAPS_AVX__ )
    typedef __m256 gaps_packed_t;
    const unsigned index_increment = 8;
    #define SET_SCALAR(x) _mm256_set1_ps(x)
    #define LOAD_PACKED(x) _mm256_load_ps(x)
    #define STORE_PACKED(p,x) _mm256_store_ps(p,x)
    #define ADD_PACKED(a,b) _mm256_add_ps(a,b)
    #define SUB_PACKED(a,b) _mm256_sub_ps(a,b)
    #define MUL_PACKED(a,b) _mm256_mul_ps(a,b)
    #define DIV_PACKED(a,b) _mm256_div_ps(a,b)
#elif defined( __GAPS_SSE__ )
    typedef __m128 gaps_packed_t;
    const unsigned index_increment = 4;
    #define SET_SCALAR(x) _mm_set1_ps(x)
    #define LOAD_PACKED(x) _mm_load_ps(x)
    #define STORE_PACKED(p,x) _mm_store_ps(p,x)
    #define ADD_PACKED(a,b) _mm_add_ps(a,b)
    #define SUB_PACKED(a,b) _mm_sub_ps(a,b)
    #define MUL_PACKED(a,b) _mm_mul_ps(a,b)
    #define DIV_PACKED(a,b) _mm_div_ps(a,b)
#else
    typedef float gaps_packed_t;
    const unsigned index_increment = 1;
    #define SET_SCALAR(x) x
    #define LOAD_PACKED(x) *(x)
    #define STORE_PACKED(p,x) *(p) = (x)
    #define ADD_PACKED(a,b) ((a)+(b))
    #define SUB_PACKED(a,b) ((a)-(b))
    #define MUL_PACKED(a,b) ((a)*(b))
    #define DIV_PACKED(a,b) ((a)/(b))
#endif

class Index
{
private:

    unsigned index;

public:

    explicit Index(unsigned i) : index(i) {}
    Index& operator=(unsigned val) { index = val; return *this; }
    bool operator<(unsigned comp) { return index < comp; }
    bool operator<=(unsigned comp) { return index <= comp; }
    void operator++() { index += index_increment; }
    unsigned value() const { return index; }
    unsigned increment() const { return index_increment; }
    friend const float* operator+(const float *ptr, Index ndx);
    friend float* operator+(float *ptr, Index ndx);
};

inline const float* operator+(const float *ptr, Index ndx) { return ptr + ndx.index; }
inline float* operator+(float *ptr, Index ndx) { return ptr + ndx.index; }

class packedFloat
{
private:

    gaps_packed_t mData;

public:

    packedFloat() : mData() {}
    explicit packedFloat(float val) : mData(SET_SCALAR(val)) {}
#if defined( __GAPS_SSE__ ) || defined( __GAPS_AVX__ ) || defined( __GAPS_AVX512__ )
    explicit packedFloat(gaps_packed_t val) : mData(val) {}
#endif

    packedFloat operator+(packedFloat b) const { return packedFloat(ADD_PACKED(mData, b.mData)); }
    packedFloat operator-(packedFloat b) const { return packedFloat(SUB_PACKED(mData, b.mData)); }
    packedFloat operator*(packedFloat b) const { return packedFloat(MUL_PACKED(mData, b.mData)); }
    packedFloat operator/(packedFloat b) const { return packedFloat(DIV_PACKED(mData, b.mData)); }

    void operator+=(packedFloat val) { mData = ADD_PACKED(mData, val.mData); }
    void load(const float *ptr) { mData = LOAD_PACKED(ptr); }
    void store(float *ptr) { STORE_PACKED(ptr, mData); }

#if defined( __GAPS_AVX__ )
    float scalar()
    {
        float* ra = reinterpret_cast<float*>(&mData); // NOLINT
        mData = _mm256_hadd_ps(mData, mData);
        mData = _mm256_hadd_ps(mData, mData);
        return ra[0] + ra[4];
    }
#elif defined( __GAPS_SSE__ )
    float scalar()
    {
        float* ra = reinterpret_cast<float*>(&mData); // NOLINT
        return ra[0] + ra[1] + ra[2] + ra[3];
    }
#else
    float scalar()
    {
        return mData;
    }
#endif
};

} // namespace simd
} // namespace gaps

#endif // __COGAPS_SIMD_H__

