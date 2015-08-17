// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cogapsMap
Rcpp::List cogapsMap(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame FixedPatt, Rcpp::DataFrame ABinsFrame, Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums);
RcppExport SEXP CoGAPS_cogapsMap(SEXP DFrameSEXP, SEXP SFrameSEXP, SEXP FixedPattSEXP, SEXP ABinsFrameSEXP, SEXP PBinsFrameSEXP, SEXP ConfigSEXP, SEXP ConfigNumsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type DFrame(DFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type SFrame(SFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type FixedPatt(FixedPattSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ABinsFrame(ABinsFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type PBinsFrame(PBinsFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type Config(ConfigSEXP );
        Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ConfigNums(ConfigNumsSEXP );
        Rcpp::List __result = cogapsMap(DFrame, SFrame, FixedPatt, ABinsFrame, PBinsFrame, Config, ConfigNums);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cogapsMapTest
Rcpp::List cogapsMapTest(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame FixedPatt, Rcpp::DataFrame ABinsFrame, Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums);
RcppExport SEXP CoGAPS_cogapsMapTest(SEXP DFrameSEXP, SEXP SFrameSEXP, SEXP FixedPattSEXP, SEXP ABinsFrameSEXP, SEXP PBinsFrameSEXP, SEXP ConfigSEXP, SEXP ConfigNumsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type DFrame(DFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type SFrame(SFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type FixedPatt(FixedPattSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ABinsFrame(ABinsFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type PBinsFrame(PBinsFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type Config(ConfigSEXP );
        Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ConfigNums(ConfigNumsSEXP );
        Rcpp::List __result = cogapsMapTest(DFrame, SFrame, FixedPatt, ABinsFrame, PBinsFrame, Config, ConfigNums);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cogaps
Rcpp::List cogaps(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame ABinsFrame, Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums);
RcppExport SEXP CoGAPS_cogaps(SEXP DFrameSEXP, SEXP SFrameSEXP, SEXP ABinsFrameSEXP, SEXP PBinsFrameSEXP, SEXP ConfigSEXP, SEXP ConfigNumsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type DFrame(DFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type SFrame(SFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ABinsFrame(ABinsFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type PBinsFrame(PBinsFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type Config(ConfigSEXP );
        Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ConfigNums(ConfigNumsSEXP );
        Rcpp::List __result = cogaps(DFrame, SFrame, ABinsFrame, PBinsFrame, Config, ConfigNums);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cogapsTest
Rcpp::List cogapsTest(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame ABinsFrame, Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums);
RcppExport SEXP CoGAPS_cogapsTest(SEXP DFrameSEXP, SEXP SFrameSEXP, SEXP ABinsFrameSEXP, SEXP PBinsFrameSEXP, SEXP ConfigSEXP, SEXP ConfigNumsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type DFrame(DFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type SFrame(SFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type ABinsFrame(ABinsFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::DataFrame >::type PBinsFrame(PBinsFrameSEXP );
        Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type Config(ConfigSEXP );
        Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ConfigNums(ConfigNumsSEXP );
        Rcpp::List __result = cogapsTest(DFrame, SFrame, ABinsFrame, PBinsFrame, Config, ConfigNums);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
