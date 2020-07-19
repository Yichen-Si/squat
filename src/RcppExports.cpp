// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// squat_single_binom_unidir
double squat_single_binom_unidir(int x, int n, double p, bool var_adj, double approx_under);
RcppExport SEXP _squat_squat_single_binom_unidir(SEXP xSEXP, SEXP nSEXP, SEXP pSEXP, SEXP var_adjSEXP, SEXP approx_underSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type var_adj(var_adjSEXP);
    Rcpp::traits::input_parameter< double >::type approx_under(approx_underSEXP);
    rcpp_result_gen = Rcpp::wrap(squat_single_binom_unidir(x, n, p, var_adj, approx_under));
    return rcpp_result_gen;
END_RCPP
}
// squat_single_binom_bidir
double squat_single_binom_bidir(int x, int n, double p, bool pos_only, bool var_adj, double approx_under);
RcppExport SEXP _squat_squat_single_binom_bidir(SEXP xSEXP, SEXP nSEXP, SEXP pSEXP, SEXP pos_onlySEXP, SEXP var_adjSEXP, SEXP approx_underSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type pos_only(pos_onlySEXP);
    Rcpp::traits::input_parameter< bool >::type var_adj(var_adjSEXP);
    Rcpp::traits::input_parameter< double >::type approx_under(approx_underSEXP);
    rcpp_result_gen = Rcpp::wrap(squat_single_binom_bidir(x, n, p, pos_only, var_adj, approx_under));
    return rcpp_result_gen;
END_RCPP
}
// squat_multi_binom_unidir
NumericVector squat_multi_binom_unidir(IntegerVector xs, NumericVector sizes, NumericVector ps, bool var_adj, double approx_under);
RcppExport SEXP _squat_squat_multi_binom_unidir(SEXP xsSEXP, SEXP sizesSEXP, SEXP psSEXP, SEXP var_adjSEXP, SEXP approx_underSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< bool >::type var_adj(var_adjSEXP);
    Rcpp::traits::input_parameter< double >::type approx_under(approx_underSEXP);
    rcpp_result_gen = Rcpp::wrap(squat_multi_binom_unidir(xs, sizes, ps, var_adj, approx_under));
    return rcpp_result_gen;
END_RCPP
}
// squat_multi_binom_bidir
NumericVector squat_multi_binom_bidir(IntegerVector xs, NumericVector sizes, NumericVector ps, bool pos_only, bool var_adj, double approx_under);
RcppExport SEXP _squat_squat_multi_binom_bidir(SEXP xsSEXP, SEXP sizesSEXP, SEXP psSEXP, SEXP pos_onlySEXP, SEXP var_adjSEXP, SEXP approx_underSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< bool >::type pos_only(pos_onlySEXP);
    Rcpp::traits::input_parameter< bool >::type var_adj(var_adjSEXP);
    Rcpp::traits::input_parameter< double >::type approx_under(approx_underSEXP);
    rcpp_result_gen = Rcpp::wrap(squat_multi_binom_bidir(xs, sizes, ps, pos_only, var_adj, approx_under));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_squat_squat_single_binom_unidir", (DL_FUNC) &_squat_squat_single_binom_unidir, 5},
    {"_squat_squat_single_binom_bidir", (DL_FUNC) &_squat_squat_single_binom_bidir, 6},
    {"_squat_squat_multi_binom_unidir", (DL_FUNC) &_squat_squat_multi_binom_unidir, 5},
    {"_squat_squat_multi_binom_bidir", (DL_FUNC) &_squat_squat_multi_binom_bidir, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_squat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}