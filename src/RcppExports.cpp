// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// min_gower_rcpp
NumericVector min_gower_rcpp(DataFrame data1_r, DataFrame data2_r);
RcppExport SEXP _flexsdm_min_gower_rcpp(SEXP data1_rSEXP, SEXP data2_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type data1_r(data1_rSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type data2_r(data2_rSEXP);
    rcpp_result_gen = Rcpp::wrap(min_gower_rcpp(data1_r, data2_r));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flexsdm_min_gower_rcpp", (DL_FUNC) &_flexsdm_min_gower_rcpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_flexsdm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
