// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cogapsMap
Rcpp::List cogapsMap(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame FixedPatt, Rcpp::DataFrame ABinsFrame, Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums, int seed, bool messages);
RcppExport SEXP _CoGAPS_cogapsMap(SEXP DFrameSEXP, SEXP SFrameSEXP, SEXP FixedPattSEXP, SEXP ABinsFrameSEXP, SEXP PBinsFrameSEXP, SEXP ConfigSEXP, SEXP ConfigNumsSEXP, SEXP seedSEXP, SEXP messagesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type DFrame(DFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type SFrame(SFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type FixedPatt(FixedPattSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ABinsFrame(ABinsFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type PBinsFrame(PBinsFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type Config(ConfigSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ConfigNums(ConfigNumsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type messages(messagesSEXP);
    rcpp_result_gen = Rcpp::wrap(cogapsMap(DFrame, SFrame, FixedPatt, ABinsFrame, PBinsFrame, Config, ConfigNums, seed, messages));
    return rcpp_result_gen;
END_RCPP
}
// cogapsMapTest
Rcpp::List cogapsMapTest(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame FixedPatt, Rcpp::DataFrame ABinsFrame, Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums, int seed);
RcppExport SEXP _CoGAPS_cogapsMapTest(SEXP DFrameSEXP, SEXP SFrameSEXP, SEXP FixedPattSEXP, SEXP ABinsFrameSEXP, SEXP PBinsFrameSEXP, SEXP ConfigSEXP, SEXP ConfigNumsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type DFrame(DFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type SFrame(SFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type FixedPatt(FixedPattSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ABinsFrame(ABinsFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type PBinsFrame(PBinsFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type Config(ConfigSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ConfigNums(ConfigNumsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(cogapsMapTest(DFrame, SFrame, FixedPatt, ABinsFrame, PBinsFrame, Config, ConfigNums, seed));
    return rcpp_result_gen;
END_RCPP
}
// cogaps
Rcpp::List cogaps(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame ABinsFrame, Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums, int seed, bool messages);
RcppExport SEXP _CoGAPS_cogaps(SEXP DFrameSEXP, SEXP SFrameSEXP, SEXP ABinsFrameSEXP, SEXP PBinsFrameSEXP, SEXP ConfigSEXP, SEXP ConfigNumsSEXP, SEXP seedSEXP, SEXP messagesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type DFrame(DFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type SFrame(SFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ABinsFrame(ABinsFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type PBinsFrame(PBinsFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type Config(ConfigSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ConfigNums(ConfigNumsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type messages(messagesSEXP);
    rcpp_result_gen = Rcpp::wrap(cogaps(DFrame, SFrame, ABinsFrame, PBinsFrame, Config, ConfigNums, seed, messages));
    return rcpp_result_gen;
END_RCPP
}
// cogapsTest
Rcpp::List cogapsTest(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame ABinsFrame, Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums, int seed);
RcppExport SEXP _CoGAPS_cogapsTest(SEXP DFrameSEXP, SEXP SFrameSEXP, SEXP ABinsFrameSEXP, SEXP PBinsFrameSEXP, SEXP ConfigSEXP, SEXP ConfigNumsSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type DFrame(DFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type SFrame(SFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ABinsFrame(ABinsFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type PBinsFrame(PBinsFrameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type Config(ConfigSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ConfigNums(ConfigNumsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(cogapsTest(DFrame, SFrame, ABinsFrame, PBinsFrame, Config, ConfigNums, seed));
    return rcpp_result_gen;
END_RCPP
}
// run_catch_unit_tests
int run_catch_unit_tests();
RcppExport SEXP _CoGAPS_run_catch_unit_tests() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(run_catch_unit_tests());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CoGAPS_cogapsMap", (DL_FUNC) &_CoGAPS_cogapsMap, 9},
    {"_CoGAPS_cogapsMapTest", (DL_FUNC) &_CoGAPS_cogapsMapTest, 8},
    {"_CoGAPS_cogaps", (DL_FUNC) &_CoGAPS_cogaps, 8},
    {"_CoGAPS_cogapsTest", (DL_FUNC) &_CoGAPS_cogapsTest, 7},
    {"_CoGAPS_run_catch_unit_tests", (DL_FUNC) &_CoGAPS_run_catch_unit_tests, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_CoGAPS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
