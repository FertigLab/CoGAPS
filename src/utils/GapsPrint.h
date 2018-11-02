#ifndef __COGAPS_GAPS_PRINT_H__
#define __COGAPS_GAPS_PRINT_H__

#ifdef __GAPS_R_BUILD__
    #define gaps_printf Rprintf
    #define gaps_cout Rcpp::Rcout
    #define gaps_flush(x) R_FlushConsole()
#else
    #define gaps_printf printf
    #define gaps_cout std::cout
    #define gaps_flush(x) fflush(stdout)
#endif

#endif // __COGAPS_GAPS_PRINT_H__
