#ifndef __COGAPS_GAPS_PRINT_H__
#define __COGAPS_GAPS_PRINT_H__

#ifdef __GAPS_R_BUILD__
    #define gaps_printf Rprintf
    #define gaps_cout Rcpp::Rcout
#else
    #define gaps_printf printf
    #define gaps_cout std::cout
#endif

#endif // __COGAPS_GAPS_PRINT_H__