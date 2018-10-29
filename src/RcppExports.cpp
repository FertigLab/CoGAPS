// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cogaps_cpp_from_file
Rcpp::List cogaps_cpp_from_file(const Rcpp::CharacterVector& data, const Rcpp::List& allParams, const Rcpp::Nullable<Rcpp::CharacterVector>& uncertainty, const Rcpp::Nullable<Rcpp::IntegerVector>& indices, const Rcpp::Nullable<Rcpp::NumericMatrix>& fixedMatrix, bool isMaster);
RcppExport SEXP _CoGAPS_cogaps_cpp_from_file(SEXP dataSEXP, SEXP allParamsSEXP, SEXP uncertaintySEXP, SEXP indicesSEXP, SEXP fixedMatrixSEXP, SEXP isMasterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type allParams(allParamsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::CharacterVector>& >::type uncertainty(uncertaintySEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::IntegerVector>& >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type fixedMatrix(fixedMatrixSEXP);
    Rcpp::traits::input_parameter< bool >::type isMaster(isMasterSEXP);
    rcpp_result_gen = Rcpp::wrap(cogaps_cpp_from_file(data, allParams, uncertainty, indices, fixedMatrix, isMaster));
    return rcpp_result_gen;
END_RCPP
}
// cogaps_cpp
Rcpp::List cogaps_cpp(const Rcpp::NumericMatrix& data, const Rcpp::List& allParams, const Rcpp::Nullable<Rcpp::NumericMatrix>& uncertainty, const Rcpp::Nullable<Rcpp::IntegerVector>& indices, const Rcpp::Nullable<Rcpp::NumericMatrix>& fixedMatrix, bool isMaster);
RcppExport SEXP _CoGAPS_cogaps_cpp(SEXP dataSEXP, SEXP allParamsSEXP, SEXP uncertaintySEXP, SEXP indicesSEXP, SEXP fixedMatrixSEXP, SEXP isMasterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type allParams(allParamsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type uncertainty(uncertaintySEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::IntegerVector>& >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type fixedMatrix(fixedMatrixSEXP);
    Rcpp::traits::input_parameter< bool >::type isMaster(isMasterSEXP);
    rcpp_result_gen = Rcpp::wrap(cogaps_cpp(data, allParams, uncertainty, indices, fixedMatrix, isMaster));
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
    {"_CoGAPS_cogaps_cpp_from_file", (DL_FUNC) &_CoGAPS_cogaps_cpp_from_file, 6},
    {"_CoGAPS_cogaps_cpp", (DL_FUNC) &_CoGAPS_cogaps_cpp, 6},
    {"_CoGAPS_getBuildReport_cpp", (DL_FUNC) &_CoGAPS_getBuildReport_cpp, 0},
    {"_CoGAPS_run_catch_unit_tests", (DL_FUNC) &_CoGAPS_run_catch_unit_tests, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_CoGAPS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
