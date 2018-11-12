#ifndef __COGAPS_SIMD_H__
#define __COGAPS_SIMD_H__

#if defined ( __AVX2__ ) || defined ( __AVX__ )

    #define SIMD_INC 8
    #define __GAPS_AVX__
    #include <immintrin.h>
    typedef __m256 gaps_packed_t;
    #define SET_SCALAR(x) _mm256_set1_ps(x)
    #define LOAD_PACKED(x) _mm256_load_ps(x)
    #define STORE_PACKED(p,x) _mm256_store_ps(p,x)
    #define ADD_PACKED(a,b) _mm256_add_ps(a,b)
    #define SUB_PACKED(a,b) _mm256_sub_ps(a,b)
    #define MUL_PACKED(a,b) _mm256_mul_ps(a,b)
    #define DIV_PACKED(a,b) _mm256_div_ps(a,b)

#elif defined ( __SSE4_2__ ) || defined ( __SSE4_1__ )

    #define SIMD_INC 4
    #define __GAPS_SSE__
    #include <nmmintrin.h>
    typedef __m128 gaps_packed_t;
    #define SET_SCALAR(x) _mm_set1_ps(x)
    #define LOAD_PACKED(x) _mm_load_ps(x)
    #define STORE_PACKED(p,x) _mm_store_ps(p,x)
    #define ADD_PACKED(a,b) _mm_add_ps(a,b)
    #define SUB_PACKED(a,b) _mm_sub_ps(a,b)
    #define MUL_PACKED(a,b) _mm_mul_ps(a,b)
    #define DIV_PACKED(a,b) _mm_div_ps(a,b)

#else

    typedef float gaps_packed_t;
    #define SIMD_INC 1
    #define SET_SCALAR(x) x
    #define LOAD_PACKED(x) *(x)
    #define STORE_PACKED(p,x) *(p) = (x)
    #define ADD_PACKED(a,b) ((a)+(b))
    #define SUB_PACKED(a,b) ((a)-(b))
    #define MUL_PACKED(a,b) ((a)*(b))
    #define DIV_PACKED(a,b) ((a)/(b))

#endif

namespace gaps
{
namespace simd
{

class Index
{
public:

    explicit Index(unsigned i) : index(i) {}
    Index& operator=(unsigned val) { index = val; return *this; }
    bool operator<(unsigned comp) const { return index < comp; }
    bool operator<=(unsigned comp) const { return index <= comp; }
    void operator++() { index += gaps::simd::Index::increment(); }
    unsigned value() const { return index; }
    
    static unsigned increment()
    {
        return SIMD_INC;
    }

    friend const float* operator+(const float *ptr, Index ndx);
    friend float* operator+(float *ptr, Index ndx);

private:

    unsigned index;
};

inline const float* operator+(const float *ptr, Index ndx) { return ptr + ndx.index; }
inline float* operator+(float *ptr, Index ndx) { return ptr + ndx.index; }

class PackedFloat
{
public:

    PackedFloat() : mData(SET_SCALAR(0.f)) {}
    explicit PackedFloat(float val) : mData(SET_SCALAR(val)) {}
#if defined( __GAPS_SSE__ ) || defined( __GAPS_AVX__ ) // avoid redefintion when gaps_packed_t == float
    explicit PackedFloat(gaps_packed_t val) : mData(val) {}
#endif

    PackedFloat operator+(PackedFloat b) const { return PackedFloat(ADD_PACKED(mData, b.mData)); }
    PackedFloat operator-(PackedFloat b) const { return PackedFloat(SUB_PACKED(mData, b.mData)); }
    PackedFloat operator*(PackedFloat b) const { return PackedFloat(MUL_PACKED(mData, b.mData)); }
    PackedFloat operator/(PackedFloat b) const { return PackedFloat(DIV_PACKED(mData, b.mData)); }

    void operator+=(PackedFloat val) { mData = ADD_PACKED(mData, val.mData); }
    void load(const float *ptr) { mData = LOAD_PACKED(ptr); }
    void store(float *ptr) { STORE_PACKED(ptr, mData); }

    float scalar()
    {
    #if defined( __GAPS_AVX__ )

        float* ra = reinterpret_cast<float*>(&mData); // NOLINT
        mData = _mm256_hadd_ps(mData, mData);
        mData = _mm256_hadd_ps(mData, mData);
        return ra[0] + ra[4];

    #elif defined( __GAPS_SSE__ )

        float* ra = reinterpret_cast<float*>(&mData); // NOLINT
        return ra[0] + ra[1] + ra[2] + ra[3];

    #else

        return mData;

    #endif
    }

private:

    gaps_packed_t mData;
};

inline float getScalar(gaps_packed_t pf)
{
    #if defined( __GAPS_AVX__ )

        pf = _mm256_hadd_ps(pf, pf);
        pf = _mm256_hadd_ps(pf, pf);
        float* ra = reinterpret_cast<float*>(&pf); // NOLINT
        return ra[0] + ra[4];

    #elif defined( __GAPS_SSE__ )

        float* ra = reinterpret_cast<float*>(&pf); // NOLINT
        return ra[0] + ra[1] + ra[2] + ra[3];

    #else

        return pf;

    #endif
}

} // namespace simd
} // namespace gaps

#endif // __COGAPS_SIMD_H__

