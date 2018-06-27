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

#if (defined(_M_AMD64) || defined(_M_X64) || defined(__amd64)) && ! defined(__x86_64__)
    #define __x86_64__ 1
#endif

#ifdef _OPENMP
    #define __GAPS_OPENMP__
#endif

// used to convert defined macro values into strings
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

inline std::string buildReport()
{
#if defined( __clang__ )
    std::string compiler = "Compiled with Clang\n";
#elif defined( __INTEL_COMPILER )
    std::string compiler = "Compiled with Intel ICC/ICPC\n";
#elif defined( __GNUC__ )
    std::string compiler = "Compiled with GCC v" + std::string(STR( __GNUC__ ))
    + "." + std::string(STR( __GNUC_MINOR__ )) + '\n';
#elif defined( _MSC_VER )
    std::string compiler = "Compiled with Microsoft Visual Studio\n";
#endif

#if defined( __GAPS_AVX__ )
    std::string simd = "AVX enabled\n";
#elif defined( __GAPS_SSE__ )
    std::string simd = "SSE enabled\n";
#else
    std::string simd = "SIMD not enabled\n";
#endif

#ifdef __GAPS_OPENMP__
    std::string openmp = "Compiled with OpenMP\n";
#else
    std::string openmp = "Compiler did not support OpenMP\n";
#endif

    return compiler + simd + openmp;
}

#endif // __COGAPS_SIMD_H__

