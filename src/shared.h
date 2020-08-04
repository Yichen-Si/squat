
#ifndef EDCPP_SHARED_H
#define EDCPP_SHARED_H

#define STRICT_R_HEADERS
#include <Rcpp.h>

using namespace Rcpp;
// Constants

#define LOGTWO 0.6931471805599452862268 // log(2.0)
#define LOGSMALL -13.81551  // log(1e-6)
#define SMALL 1e-6
#define MIN_SD 0.01 // minimum SD

static const double SQRT_2_PI    = 2.506628274631000241612;  // sqrt(2*pi)
static const double PHI_0        = 0.3989422804014327028632; // dnorm(0)
static const double LOG_2F       = 0.6931471805599452862268; // log(2)

static const double MIN_DIFF_EPS = 1e-8;

// MACROS

#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))
#define LOGDIFF(a,b) ((a > b) ? a + log(1.0-exp(b-a)) : b + log(1.0-exp(a-b))) // |a-b| in logscale
#define LOGADD(a,b) ((a < b) ? b + log(1.0+exp(a-b)) : a + log(1.0+exp(b-a))   // a+b in logscale

// functions

bool isInteger(double x, bool warn = true);
double finite_max_int(const Rcpp::NumericVector& x);
double rng_unif();         // standard uniform

NumericVector cpp_dbbinom(const NumericVector& x, const NumericVector& size,
    const NumericVector& alpha, const NumericVector& beta,
    const bool& log_prob = false);
NumericVector cpp_pbbinom(const NumericVector& x, const NumericVector& size,
    const NumericVector& alpha, const NumericVector& beta,
    const bool& lower_tail = true, const bool& log_prob = false);
double cpp_pbbinom_one(double x, double size, double alpha, double beta,
    bool lower_tail = true, bool log_prob = false);
double cpp_dbbinom_one(double x, double size, double alpha, double beta,
                       bool log_prob = false);

// inline functions

inline bool tol_equal(double x, double y);
inline double phi(double x);
inline double lphi(double x);
inline double Phi(double x);
inline double InvPhi(double x);
inline double factorial(double x);
inline double lfactorial(double x);
inline double rng_sign();
inline bool is_large_int(double x);
inline double to_dbl(int x);
inline int to_pos_int(double x);
inline double trunc_p(double x);

#include "shared_inline.h"


#endif
