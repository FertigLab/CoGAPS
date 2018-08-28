#ifndef __COGAPS_GAPS_ASSERT_H__
#define __COGAPS_GAPS_ASSERT_H__

#include "GapsPrint.h"

#ifdef __GAPS_R_BUILD__
#include <Rcpp.h>
#endif

#ifndef __GAPS_R_BUILD__
#include <cstdlib>
#include <iostream>
#endif

#ifdef __GAPS_R_BUILD__
#define gaps_stop() Rcpp::stop("CoGAPS terminated")
#else
#define gaps_stop() std::exit(0)
#endif

// NOLINTNEXTLINE
#define GAPS_ERROR(msg) do {gaps_cout << "error: " << msg << '\n'; gaps_stop();} while(0)

#ifdef GAPS_DEBUG
    #define GAPS_ASSERT(cond)                                             \
        do {                                                              \
            if (!(cond))                                                  \
            {                                                             \
                gaps_printf("assert failed %s %d\n", __FILE__, __LINE__); \
                gaps_stop();                                              \
            }                                                             \
        } while(0)

    #define GAPS_ASSERT_MSG(cond, msg)                                  \
        do {                                                            \
            if (!(cond))                                                \
            {                                                           \
                gaps_cout << "assert failed " << __FILE__ << " " <<     \
                    __LINE__ << ", " << msg << '\n';                   \
                gaps_stop();                                            \
            }                                                           \
        } while(0)

    #define DEBUG_PING gaps_printf("here %s %d\n", __FILE__, __LINE__);
#else
    #define GAPS_ASSERT(cond) do {} while(0)
    #define GAPS_ASSERT_MSG(cond, msg) do {} while(0)
    #define DEBUG_PING   
#endif

#endif // __COGAPS_GAPS_ASSERT_H__