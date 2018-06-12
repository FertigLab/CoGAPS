#include "GapsDispatcher.h"

#include <string>

Rcpp::List GapsDispatcher::run()
{
    for (unsigned i = 0; i < mMaxIterations; ++i)
    {
        // run all GapsRunners for K iterations

        // sync P matrices

        // broadcast master P and weight
    }

    // sync all P matrices to same master

    for (unsigned i = 0; i < mMaxIterations; ++i)
    {
        // runner all samplers with sampling turned on
    }

    // stitch together matrices
}