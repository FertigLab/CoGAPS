#ifndef __COGAPS_GAPS_ASSERT_H__
#define __COGAPS_GAPS_ASSERT_H__

#include "GapsPrint.h"

#define GAPS_REAL_ASSERT

#ifdef GAPS_DEBUG
#define GAPS_REAL_ASSERT
#endif

#ifdef __GAPS_R_BUILD__
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Rcpp.h>
//#pragma GCC diagnostic pop
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

#ifdef GAPS_REAL_ASSERT
    
    #define GAPS_ERROR(msg)                                      \
        {                                                        \
            std::cout << "error: " << msg << '\n'                \
                << __FILE__ << __LINE__ << std::flush;           \
            gaps_stop();                                         \
        }


    #define GAPS_ASSERT(cond)                                     \
        if (!(cond))                                              \
        {                                                         \
            std::cout<< "GAPS assert failed \nat " <<         \
            __FILE__":" << __LINE__ << '\n' << std::flush;             \
            gaps_stop();                                          \
        }                                                         \

    #define GAPS_ASSERT_MSG(cond, msg)                              \
        if (!(cond))                                                \
        {                                                           \
            std::cout << msg <<"\nat "<<__FILE__ << ":" <<     \
                __LINE__ << '\n' << std::flush;      \
            gaps_stop();                                            \
        }                                                          

    #define DEBUG_PING gaps_printf("here %s %d\n", __FILE__, __LINE__);
#else
    #define GAPS_ASSERT(cond) do {} while(0)
    #define GAPS_ASSERT_MSG(cond, msg) do {} while(0)
    #define DEBUG_PING   
#endif

#endif // __COGAPS_GAPS_ASSERT_H__
