// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cogaps_from_file_cpp
Rcpp::List cogaps_from_file_cpp(const Rcpp::CharacterVector& data, const Rcpp::List& allParams, const Rcpp::Nullable<Rcpp::CharacterVector>& uncertainty);
RcppExport SEXP _CoGAPS_cogaps_from_file_cpp(SEXP dataSEXP, SEXP allParamsSEXP, SEXP uncertaintySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type allParams(allParamsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::CharacterVector>& >::type uncertainty(uncertaintySEXP);
    rcpp_result_gen = Rcpp::wrap(cogaps_from_file_cpp(data, allParams, uncertainty));
    return rcpp_result_gen;
END_RCPP
}
// cogaps_cpp
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix& data, const Rcpp::List& allParams, const Rcpp::Nullable<Rcpp::NumericMatrix>& uncertainty);
RcppExport SEXP _CoGAPS_cogaps_cpp(SEXP dataSEXP, SEXP allParamsSEXP, SEXP uncertaintySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type allParams(allParamsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type uncertainty(uncertaintySEXP);
    rcpp_result_gen = Rcpp::wrap(cogaps_cpp(data, allParams, uncertainty));
    return rcpp_result_gen;
END_RCPP
}
// getBuildReport_cpp
std::string getBuildReport_cpp();
RcppExport SEXP _CoGAPS_getBuildReport_cpp() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getBuildReport_cpp());
    return rcpp_result_gen;
END_RCPP
}
// checkpointsEnabled_cpp
bool checkpointsEnabled_cpp();
RcppExport SEXP _CoGAPS_checkpointsEnabled_cpp() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(checkpointsEnabled_cpp());
    return rcpp_result_gen;
END_RCPP
}
// compiledWithOpenMPSupport_cpp
bool compiledWithOpenMPSupport_cpp();
RcppExport SEXP _CoGAPS_compiledWithOpenMPSupport_cpp() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(compiledWithOpenMPSupport_cpp());
    return rcpp_result_gen;
END_RCPP
}
// getFileInfo_cpp
Rcpp::List getFileInfo_cpp(const std::string& path);
RcppExport SEXP _CoGAPS_getFileInfo_cpp(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(getFileInfo_cpp(path));
    return rcpp_result_gen;
END_RCPP
}
// run_catch_unit_tests
int run_catch_unit_tests(Rcpp::String reporter);
RcppExport SEXP _CoGAPS_run_catch_unit_tests(SEXP reporterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type reporter(reporterSEXP);
    rcpp_result_gen = Rcpp::wrap(run_catch_unit_tests(reporter));
    return rcpp_result_gen;
END_RCPP
}
// run_catch_unit_tests_by_tag
int run_catch_unit_tests_by_tag(Rcpp::String tag, Rcpp::String reporter);
RcppExport SEXP _CoGAPS_run_catch_unit_tests_by_tag(SEXP tagSEXP, SEXP reporterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type tag(tagSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type reporter(reporterSEXP);
    rcpp_result_gen = Rcpp::wrap(run_catch_unit_tests_by_tag(tag, reporter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CoGAPS_cogaps_from_file_cpp", (DL_FUNC) &_CoGAPS_cogaps_from_file_cpp, 3},
    {"_CoGAPS_cogaps_cpp", (DL_FUNC) &_CoGAPS_cogaps_cpp, 3},
    {"_CoGAPS_getBuildReport_cpp", (DL_FUNC) &_CoGAPS_getBuildReport_cpp, 0},
    {"_CoGAPS_checkpointsEnabled_cpp", (DL_FUNC) &_CoGAPS_checkpointsEnabled_cpp, 0},
    {"_CoGAPS_compiledWithOpenMPSupport_cpp", (DL_FUNC) &_CoGAPS_compiledWithOpenMPSupport_cpp, 0},
    {"_CoGAPS_getFileInfo_cpp", (DL_FUNC) &_CoGAPS_getFileInfo_cpp, 1},
    {"_CoGAPS_run_catch_unit_tests", (DL_FUNC) &_CoGAPS_run_catch_unit_tests, 1},
    {"_CoGAPS_run_catch_unit_tests_by_tag", (DL_FUNC) &_CoGAPS_run_catch_unit_tests_by_tag, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_CoGAPS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
