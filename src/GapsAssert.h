#ifndef __GAPS_ASSERT_H__
#define __GAPS_ASSERT_H__

#ifdef GAPS_DEBUG
    #define GAPS_ASSERT(cond)                                           \
        do {                                                            \
            if (!(cond))                                                \
            {                                                           \
                Rprintf("assert failed %s %d\n", __FILE__, __LINE__);   \
                std::exit(0);                                           \
            }                                                           \
        } while(0)
#else
    #define GAPS_ASSERT(cond) ((void)sizeof(cond))
#endif 

#endif