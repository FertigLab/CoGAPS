#ifndef __COGAPS_GLOBAL_CONFIG_H__
#define __COGAPS_GLOBAL_CONFIG_H__

#include "../math/SIMD.h"

#include <string>

#if (defined(_M_AMD64) || defined(_M_X64) || defined(__amd64)) && ! defined(__x86_64__)
    #define __x86_64__ 1
#endif

#ifdef _OPENMP
    #define __GAPS_OPENMP__
#endif

#ifdef __GAPS_R_BUILD__
    #define gaps_check_interrupt(x) Rcpp::checkUserInterrupt(x)
#else
    #define gaps_check_interrupt(x) do {} while(0)
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
    std::string simd = "SIMD: AVX instructions enabled\n";
#elif defined( __GAPS_SSE__ )
    std::string simd = "SIMD: SSE instructions enabled\n";
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

#endif // __COGAPS_GLOBAL_CONFIG_H__