// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// f
double f(double x);
RcppExport SEXP _SC19090_f(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(f(x));
    return rcpp_result_gen;
END_RCPP
}
// randomMHC
List randomMHC(double sigma, double x0, int N);
RcppExport SEXP _SC19090_randomMHC(SEXP sigmaSEXP, SEXP x0SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(randomMHC(sigma, x0, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SC19090_f", (DL_FUNC) &_SC19090_f, 1},
    {"_SC19090_randomMHC", (DL_FUNC) &_SC19090_randomMHC, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SC19090(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}