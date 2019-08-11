// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// binaryCodingCpp
NumericMatrix binaryCodingCpp(DataFrame data);
RcppExport SEXP _Cookie_binaryCodingCpp(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(binaryCodingCpp(data));
    return rcpp_result_gen;
END_RCPP
}
// hammingCodingCpp
NumericMatrix hammingCodingCpp(DataFrame data);
RcppExport SEXP _Cookie_hammingCodingCpp(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(hammingCodingCpp(data));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Cookie_binaryCodingCpp", (DL_FUNC) &_Cookie_binaryCodingCpp, 1},
    {"_Cookie_hammingCodingCpp", (DL_FUNC) &_Cookie_hammingCodingCpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_Cookie(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
