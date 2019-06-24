#ifndef __COGAPS_GAPS_PRINT_H__
#define __COGAPS_GAPS_PRINT_H__

// printing with Rcpp doesn't display when executing inside BiocParallel

//#ifdef __GAPS_R_BUILD__

    //#include <Rcpp.h>

    //#define gaps_printf Rprintf
    //#define gaps_cout Rcpp::Rcout
    //#define gaps_flush(dummy) R_FlushConsole()

//#else

    #include <cstdio>
    #include <iostream>

    #define gaps_printf printf
    #define gaps_cout std::cout
    #define gaps_flush(dummy) fflush(stdout)

//#endif

#endif // __COGAPS_GAPS_PRINT_H__
