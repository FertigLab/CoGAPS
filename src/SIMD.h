#ifndef __COGAPS_SIMD_H__
#define __COGAPS_SIMD_H__

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#pragma message "GCC Version: " STR(__GNUC__) "." STR(__GNUC_MINOR__)

#if (defined(_M_AMD64) || defined(_M_X64) || defined(__amd64)) && ! defined(__x86_64__)
    #define __x86_64__ 1
#endif

#ifndef SSE_INSTR_SET
    #ifndef SIMD
        #define SSE_INSTR_SET 0
    #elif defined ( __AVX2__ )
        #define SSE_INSTR_SET 8
    #elif defined ( __AVX__ )
        #define SSE_INSTR_SET 7
    #elif defined ( __SSE4_2__ )
        #define SSE_INSTR_SET 6
    #elif defined ( __SSE4_1__ )
        #define SSE_INSTR_SET 5
    #elif defined ( __SSSE3__ )
        #define SSE_INSTR_SET 4
    #elif defined ( __SSE3__ )
        #define SSE_INSTR_SET 3
    #elif defined ( __SSE2__ ) || defined ( __x86_64__ )
        #define SSE_INSTR_SET 2
    #elif defined ( __SSE__ )
        #define SSE_INSTR_SET 1
    #elif defined ( _M_IX86_FP )  // Defined in MS compiler on 32bits system. 1: SSE, 2: SSE2
        #define SSE_INSTR_SET _M_IX86_FP
    #else
        #error "SIMD not supported"
        #define SSE_INSTR_SET 0
    #endif // instruction set defines
#endif // SSE_INSTR_SET

// Include the appropriate header file for intrinsic functions
#if SSE_INSTR_SET > 7                  // AVX2 and later
    #pragma message "AVX2 enabled"
    #ifdef __GNUC__
        #include <x86intrin.h>         // x86intrin.h includes header files for whatever instruction
                                       // sets are specified on the compiler command line, such as:
                                       // xopintrin.h, fma4intrin.h
    #else
        #include <immintrin.h>         // MS version of immintrin.h covers AVX, AVX2 and FMA3
    #endif // __GNUC__
#elif SSE_INSTR_SET == 7
    #pragma message "AVX enabled"
    #include <immintrin.h>
#elif SSE_INSTR_SET == 6
    #pragma message "SSE4.2 enabled"
    #include <nmmintrin.h>
#elif SSE_INSTR_SET == 5
    #pragma message "SSE4.1 enabled"
    #include <smmintrin.h>
#elif SSE_INSTR_SET == 4
    #pragma message "SSSE3 enabled"
    #include <tmmintrin.h>
#elif SSE_INSTR_SET == 3
    #pragma message "SSE3 enabled"
    #include <pmmintrin.h>
#elif SSE_INSTR_SET == 2
    #pragma message "SSE2 enabled"
    #include <emmintrin.h>
#elif SSE_INSTR_SET == 1
    #pragma message "SSE enabled"
    #include <xmmintrin.h>
#endif

namespace gaps
{
namespace simd
{

#if SSE_INSTR_SET == 0
    const unsigned index_increment = 1;
#elif SSE_INSTR_SET == 7
    const unsigned index_increment = 8;
#else
    const unsigned index_increment = 1;
#endif

class Index
{
private:
    unsigned index;
public:
    Index(unsigned i) : index(i) {}
    void operator=(unsigned val) { index = val; }
    bool operator<(unsigned comp) { return index < comp; }
    void operator++() { index += index_increment; }
    friend const float* operator+(const float *ptr, Index ndx);
};

inline const float* operator+(const float *ptr, Index ndx)
{
    return ptr + ndx.index;
}

#ifdef SIMD

#include <immintrin.h>

class packedFloat
{
private:
    __m128 mVec;
public:
    vec4f() {}
    vec4f(float val) : mVec(_mm_load_ps1(&val)) {}
    vec4f(const float *ptr) : mVec(_mm_load_ps(ptr)) {}
    vec4f(__m128 val) : mVec(val) {}
    void operator=(float val) { mVec = _mm_load_ps1(&val); }
    void operator+=(vec4f val) { mVec = mVec + val.mVec; }

    void load(const float* ptr) { mVec = _mm_load_ps(ptr); }

    vec4f operator+(vec4f b) const { return _mm_add_ps(mVec, b.mVec); }
    vec4f operator-(vec4f b) const { return _mm_sub_ps(mVec, b.mVec); }
    vec4f operator*(vec4f b) const { return _mm_mul_ps(mVec, b.mVec); }
    vec4f operator/(vec4f b) const { return _mm_div_ps(mVec, b.mVec); }

    friend float scalar(vec4f val);
};

inline float scalar(vec4f val)
{
    val.mVec = _mm_hadd_ps(val.mVec, val.mVec);
    val.mVec = _mm_hadd_ps(val.mVec, val.mVec);
    float ret;
    _mm_store_ss(&ret, val.mVec);
    return ret;
}

#else
#warning "SIMD not being used"

class packedFloat
{
private:
    float mVec;
public:
    vec4f() {}
    vec4f(float val) : mVec(val) {}
    vec4f(const float *ptr) : mVec(*ptr) {}
    void operator=(float val) { mVec = val; }
    void operator+=(vec4f val) { mVec += val.mVec; }

    vec4f operator+(vec4f b) { return mVec + b.mVec; }
    vec4f operator-(vec4f b) { return mVec - b.mVec; }
    vec4f operator*(vec4f b) { return mVec * b.mVec; }
    vec4f operator/(vec4f b) { return mVec / b.mVec; }

    friend float scalar(vec4f val);
};

inline vec4f load(const float* ptr) { return *ptr; }
inline float scalar(vec4f val) { return val.mVec; }

#endif // SIMD

} // namespace simd
} // namespace gaps

#endif // __COGAPS_SIMD_H__
