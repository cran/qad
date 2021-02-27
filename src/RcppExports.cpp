// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// build_checkerboard_weights
NumericMatrix build_checkerboard_weights(const NumericVector& X, const NumericVector& Y, R_xlen_t resolution);
RcppExport SEXP _qad_build_checkerboard_weights(SEXP XSEXP, SEXP YSEXP, SEXP resolutionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type resolution(resolutionSEXP);
    rcpp_result_gen = Rcpp::wrap(build_checkerboard_weights(X, Y, resolution));
    return rcpp_result_gen;
END_RCPP
}
// local_kernel_integral
double local_kernel_integral(const NumericMatrix& A, R_xlen_t x, R_xlen_t y, R_xlen_t N, double y_sum);
RcppExport SEXP _qad_local_kernel_integral(SEXP ASEXP, SEXP xSEXP, SEXP ySEXP, SEXP NSEXP, SEXP y_sumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A(ASEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type x(xSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type y(ySEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type y_sum(y_sumSEXP);
    rcpp_result_gen = Rcpp::wrap(local_kernel_integral(A, x, y, N, y_sum));
    return rcpp_result_gen;
END_RCPP
}
// D1_Pi
double D1_Pi(const NumericMatrix& A, R_xlen_t resolution);
RcppExport SEXP _qad_D1_Pi(SEXP ASEXP, SEXP resolutionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A(ASEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type resolution(resolutionSEXP);
    rcpp_result_gen = Rcpp::wrap(D1_Pi(A, resolution));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_qad_build_checkerboard_weights", (DL_FUNC) &_qad_build_checkerboard_weights, 3},
    {"_qad_local_kernel_integral", (DL_FUNC) &_qad_local_kernel_integral, 5},
    {"_qad_D1_Pi", (DL_FUNC) &_qad_D1_Pi, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_qad(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
