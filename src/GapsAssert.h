#ifndef __COGAPS_GAPS_ASSERT_H__
#define __COGAPS_GAPS_ASSERT_H__

#ifdef GAPS_DEBUG
    #define GAPS_ASSERT(cond)                                           \
        do {                                                            \
            if (!(cond))                                                \
            {                                                           \
                Rprintf("assert failed %s %d\n", __FILE__, __LINE__);   \
                std::exit(0);                                           \
            }                                                           \
        } while(0)

    #define GAPS_ASSERT_MSG(cond, msg)                                  \
        do {                                                            \
            if (!(cond))                                                \
            {                                                           \
                Rcpp::Rcout << "assert failed " << __FILE__ << " " <<   \
                    __LINE__ << " , " << msg << '\n';                   \
                std::exit(0);                                           \
            }                                                           \
        } while(0)
#else
    #define GAPS_ASSERT(cond) ((void)sizeof(cond))
    #define GAPS_ASSERT_MSG(cond, msg) ((void)sizeof(cond))
#endif 

#endif