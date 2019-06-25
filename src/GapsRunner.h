#ifndef __COGAPS_GAPS_RUNNER_H__
#define __COGAPS_GAPS_RUNNER_H__

struct GapsResult;
struct GapsParameters;
class Matrix;
class GapsRandomState;

#include <string>

// these two functions are the top-level functions exposed to the C++
// code that is being wrapped by any given language

namespace gaps
{
    // data stored in matrix
    GapsResult run(const Matrix &data, GapsParameters &params,
        const Matrix &uncertainty, GapsRandomState *randState);

    // data stored in file
    GapsResult run(const std::string &data, GapsParameters &params,
        const std::string &uncertainty, GapsRandomState *randState);

}; // namespace gaps

#endif // __COGAPS_GAPS_RUNNER_H__