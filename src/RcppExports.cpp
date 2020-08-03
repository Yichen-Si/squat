// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// squat_single_binom_unidir
NumericVector squat_single_binom_unidir(int x, int n, double p, bool var_adj, double approx_under);
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
NumericVector squat_single_binom_bidir(int x, int n, double p, bool pos_only, bool var_adj, double approx_under);
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
DataFrame squat_multi_binom_unidir(IntegerVector xs, NumericVector sizes, NumericVector ps, bool var_adj, double approx_under);
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
DataFrame squat_multi_binom_bidir(IntegerVector xs, NumericVector sizes, NumericVector ps, bool pos_only, bool var_adj, double approx_under);
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
// expt_truncated_normal_from_qt
double expt_truncated_normal_from_qt(double a, double b, double mu, double sd, bool lg, bool lower);
RcppExport SEXP _squat_expt_truncated_normal_from_qt(SEXP aSEXP, SEXP bSEXP, SEXP muSEXP, SEXP sdSEXP, SEXP lgSEXP, SEXP lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type lg(lgSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(expt_truncated_normal_from_qt(a, b, mu, sd, lg, lower));
    return rcpp_result_gen;
END_RCPP
}
// expt_truncated_bidir_normal_from_qt
double expt_truncated_bidir_normal_from_qt(double ql, double qu, double mu, double sd, bool lg, bool lower);
RcppExport SEXP _squat_expt_truncated_bidir_normal_from_qt(SEXP qlSEXP, SEXP quSEXP, SEXP muSEXP, SEXP sdSEXP, SEXP lgSEXP, SEXP lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< double >::type qu(quSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< bool >::type lg(lgSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(expt_truncated_bidir_normal_from_qt(ql, qu, mu, sd, lg, lower));
    return rcpp_result_gen;
END_RCPP
}
// multi_bidir_zs_from_qt_n
NumericVector multi_bidir_zs_from_qt_n(NumericVector ql, NumericVector qu, bool lg, bool lower);
RcppExport SEXP _squat_multi_bidir_zs_from_qt_n(SEXP qlSEXP, SEXP quSEXP, SEXP lgSEXP, SEXP lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type qu(quSEXP);
    Rcpp::traits::input_parameter< bool >::type lg(lgSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(multi_bidir_zs_from_qt_n(ql, qu, lg, lower));
    return rcpp_result_gen;
END_RCPP
}
// bidir_etn_var_from_qt
double bidir_etn_var_from_qt(NumericVector qt, bool lg, bool lower);
RcppExport SEXP _squat_bidir_etn_var_from_qt(SEXP qtSEXP, SEXP lgSEXP, SEXP lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type qt(qtSEXP);
    Rcpp::traits::input_parameter< bool >::type lg(lgSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(bidir_etn_var_from_qt(qt, lg, lower));
    return rcpp_result_gen;
END_RCPP
}
// binom_multi_bidir_n
List binom_multi_bidir_n(IntegerVector xs, NumericVector sizes, NumericVector ps, bool pos_only, bool var_adj, double approx_under);
RcppExport SEXP _squat_binom_multi_bidir_n(SEXP xsSEXP, SEXP sizesSEXP, SEXP psSEXP, SEXP pos_onlySEXP, SEXP var_adjSEXP, SEXP approx_underSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< bool >::type pos_only(pos_onlySEXP);
    Rcpp::traits::input_parameter< bool >::type var_adj(var_adjSEXP);
    Rcpp::traits::input_parameter< double >::type approx_under(approx_underSEXP);
    rcpp_result_gen = Rcpp::wrap(binom_multi_bidir_n(xs, sizes, ps, pos_only, var_adj, approx_under));
    return rcpp_result_gen;
END_RCPP
}
// expt_truncated_gamma_from_qt
double expt_truncated_gamma_from_qt(double a, double b, double alpha, double beta, bool lg, bool lower);
RcppExport SEXP _squat_expt_truncated_gamma_from_qt(SEXP aSEXP, SEXP bSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP lgSEXP, SEXP lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type lg(lgSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(expt_truncated_gamma_from_qt(a, b, alpha, beta, lg, lower));
    return rcpp_result_gen;
END_RCPP
}
// expt_truncated_bidir_gamma_from_qt
double expt_truncated_bidir_gamma_from_qt(double ql, double qu, double alpha, double beta, bool lg, bool lower);
RcppExport SEXP _squat_expt_truncated_bidir_gamma_from_qt(SEXP qlSEXP, SEXP quSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP lgSEXP, SEXP lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< double >::type qu(quSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type lg(lgSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(expt_truncated_bidir_gamma_from_qt(ql, qu, alpha, beta, lg, lower));
    return rcpp_result_gen;
END_RCPP
}
// squat_single_binom_unidir_g
NumericVector squat_single_binom_unidir_g(int x, int n, double p, double alpha, double beta, bool var_adj, double approx_under, bool lower);
RcppExport SEXP _squat_squat_single_binom_unidir_g(SEXP xSEXP, SEXP nSEXP, SEXP pSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP var_adjSEXP, SEXP approx_underSEXP, SEXP lowerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type var_adj(var_adjSEXP);
    Rcpp::traits::input_parameter< double >::type approx_under(approx_underSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    rcpp_result_gen = Rcpp::wrap(squat_single_binom_unidir_g(x, n, p, alpha, beta, var_adj, approx_under, lower));
    return rcpp_result_gen;
END_RCPP
}
// squat_multi_binom_dir_g
List squat_multi_binom_dir_g(IntegerVector xs, NumericVector sizes, NumericVector ps, NumericVector ys, bool var_adj, double approx_under);
RcppExport SEXP _squat_squat_multi_binom_dir_g(SEXP xsSEXP, SEXP sizesSEXP, SEXP psSEXP, SEXP ysSEXP, SEXP var_adjSEXP, SEXP approx_underSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ps(psSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< bool >::type var_adj(var_adjSEXP);
    Rcpp::traits::input_parameter< double >::type approx_under(approx_underSEXP);
    rcpp_result_gen = Rcpp::wrap(squat_multi_binom_dir_g(xs, sizes, ps, ys, var_adj, approx_under));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_squat_squat_single_binom_unidir", (DL_FUNC) &_squat_squat_single_binom_unidir, 5},
    {"_squat_squat_single_binom_bidir", (DL_FUNC) &_squat_squat_single_binom_bidir, 6},
    {"_squat_squat_multi_binom_unidir", (DL_FUNC) &_squat_squat_multi_binom_unidir, 5},
    {"_squat_squat_multi_binom_bidir", (DL_FUNC) &_squat_squat_multi_binom_bidir, 6},
    {"_squat_expt_truncated_normal_from_qt", (DL_FUNC) &_squat_expt_truncated_normal_from_qt, 6},
    {"_squat_expt_truncated_bidir_normal_from_qt", (DL_FUNC) &_squat_expt_truncated_bidir_normal_from_qt, 6},
    {"_squat_multi_bidir_zs_from_qt_n", (DL_FUNC) &_squat_multi_bidir_zs_from_qt_n, 4},
    {"_squat_bidir_etn_var_from_qt", (DL_FUNC) &_squat_bidir_etn_var_from_qt, 3},
    {"_squat_binom_multi_bidir_n", (DL_FUNC) &_squat_binom_multi_bidir_n, 6},
    {"_squat_expt_truncated_gamma_from_qt", (DL_FUNC) &_squat_expt_truncated_gamma_from_qt, 6},
    {"_squat_expt_truncated_bidir_gamma_from_qt", (DL_FUNC) &_squat_expt_truncated_bidir_gamma_from_qt, 6},
    {"_squat_squat_single_binom_unidir_g", (DL_FUNC) &_squat_squat_single_binom_unidir_g, 8},
    {"_squat_squat_multi_binom_dir_g", (DL_FUNC) &_squat_squat_multi_binom_dir_g, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_squat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
