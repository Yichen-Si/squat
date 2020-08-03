#include <Rcpp.h>
//#include <cmath>
using namespace Rcpp;

#define LOGDIFF(a,b) ((a > b) ? a + log(1.0-exp(b-a)) : b + log(1.0-exp(a-b))) // |a-b| in logscale
#define LOGADD(a,b) ((a < b) ? b + log(1.0+exp(a-b)) : a + log(1.0+exp(b-a))   // a+b in logscale
#define LOGTWO 0.6931471805599452862268 // log(2.0)
#define LOGSMALL -13.81551  // log(1e-6)
#define SMALL 1e-6;
#define MIN_SD 0.01 // minimum SD

//' @title expt_truncated_normal_from_qt
//' @description
//' Calculate expectation from truncated Normal
//'
//' @param a lower quantile
//' @param b upper quantile
//' @param mu mean of the Normal
//' @param sd standard deviation of the Normal
//' @param lg compute log-scale
//' @param lower input quantiles are lower tail
//' @return A single expected z-score for directional test
// [[Rcpp::export]]
double expt_truncated_normal_from_qt(double a, double b, 
                                     double mu=0, double sd=1, 
                                     bool lg=0, bool lower=1) {
  double denom;
  if (lg) {
    a = std::max(LOGSMALL, a); b = std::max(LOGSMALL, b);
    a = std::min(SMALL, a); b = std::min(SMALL, b);
  } else {
    a = std::min(1-SMALL,a); b = std::min(1-SMALL,b);
    a = std::max(SMALL,a); b = std::max(SMALL,b);
  }
  if (!lower) {
    if (lg) {
      b = LOGDIFF(b,0);
      a = LOGDIFF(a,0);
    } else {
      b = 1-b;
      a = 1-a;
    }
    lower = 1;
  }
  if (lg) {
    if ( LOGDIFF(a,b) < LOGSMALL ) {return R::qnorm(a,mu,sd,lower,lg);}
    denom = ( (b>a) ? exp(LOGDIFF(a,b)) : -exp(LOGDIFF(a,b)) );
  } else {
    if ( abs(a-b) < 1e-6 ) {
      return R::qnorm(a,mu,sd,lower,lg);
      }
    denom = b-a;
  }
  return mu + sd * (R::dnorm(R::qnorm(a,mu,sd,1,lg),mu,sd,0) -
                    R::dnorm(R::qnorm(b,mu,sd,1,lg),mu,sd,0) ) / denom;
}


//' @title expt_truncated_bidir_normal_from_qt
//' @description
//' Calculate expectation from truncated Normal - folded
//'
//' @param ql lower quantile
//' @param qu upper quantile
//' @param mu mean of the Normal
//' @param sd standard deviation of the Normal
//' @param lg compute log-scale
//' @param lower input quantiles are lower tail
//' @return A single expected z-score for overdispersion test
// [[Rcpp::export]]
double expt_truncated_bidir_normal_from_qt( double ql, double qu, 
                                            double mu=0, double sd=1, 
                                            bool lg=0, bool lower=1 ) {
  double ez;
  if (lg) {
    double ql2 = log(2.)+ql;
    double qu2 = log(2.)+qu;
    if (!lower) {
      double tmp = LOGDIFF(0,qu);
      qu = LOGDIFF(0,ql);
      ql = tmp;
      lower=1;
    }
    if (ql > log(0.5)) {
      ez = expt_truncated_normal_from_qt( LOGDIFF(qu2,log(2.)), LOGDIFF(ql2,log(2.)), mu,sd, lg, lower );
    } else if (qu <= log(0.5)) {
      ez = expt_truncated_normal_from_qt(ql2, qu2, mu,sd, lg, lower );
    } else {
      double w = ( 0.5-exp(ql) ) / ( exp(qu)-exp(ql) );
      ez = w * expt_truncated_normal_from_qt( ql2, 0, mu,sd, lg, lower ) +
        (1.-w) * expt_truncated_normal_from_qt( LOGDIFF(qu2,log(2.)), 0, mu,sd, lg, lower );
    }
  } else {
    if (!lower) {
      double tmp = 1-qu;
      qu = 1-ql;
      ql = tmp;
      lower=1;
    }
    if (ql > 0.5) {
      ez = expt_truncated_normal_from_qt( 2.-2.*qu, 2.-2.*ql, mu,sd, lg, lower );
    } else if (qu <= 0.5) {
      ez = expt_truncated_normal_from_qt( 2.*ql, 2.*qu, mu,sd, lg, lower );
    } else {
      double w = ( 0.5-ql ) / (qu - ql);
      ez = w * expt_truncated_normal_from_qt(2.*ql, 1, mu,sd, lg, lower ) +
        (1.-w) * expt_truncated_normal_from_qt(2.-2.*qu, 1, mu,sd, lg, lower );
    }
  }
  return ez;
}


//' @title multi_bidir_zs_from_qt_n
//' @description
//' A function to generate bi-directional z scores based on input quantiles (from an arbitrary distribution)
//'
//' @param ql A vector of lower quantiles
//' @param ql A vector of upper quantiles
//' @param lg compute log-scale
//' @param lower input quantiles are lower tail
//' @return A vector including z-scores for the input quantile pairs
// [[Rcpp::export]]
NumericVector multi_bidir_zs_from_qt_n(NumericVector ql, NumericVector qu, bool lg=0, bool lower=1) {
  int n = ql.size();
  if ( qu.size() != n ) {
    stop( "sizes of lower and upper quantiles must have the same length." );
  }
  NumericVector zs(n);
  for (int i = 0; i < n; i++) {
    zs[i] = expt_truncated_bidir_normal_from_qt( ql[i], qu[i], 0, 1, lg, lower);
  }
  return zs;
}


//' @title bidir_etn_var_from_qt
//' @description
//' A function to generate (approximated) variance for expectation based bi-directional z scores 
//' given an arbitrary distribution with the majority of mass captured by the input quantiles
//' 
//' @param qt A vector of quantiles
//' @param lg compute log-scale
//' @param lower input quantiles are lower tail
//' @return A vector including z-scores for the input quantile pairs
// [[Rcpp::export]]
double bidir_etn_var_from_qt(NumericVector qt, bool lg=0, bool lower=1) {
  int n = qt.size();
  double vz = 0.0;
  if (n > 1) {
    if (qt[n-1] < 0.5 || qt[0] > 0.5) {
      stop("The exact calculation interval must cover 0.5.");
    }
    double z;
    for (int i = 1; i < n; i++) {
      double den = lg ? exp(LOGDIFF(qt[i-1],qt[i])) : abs(qt[i-1]-qt[i]);
      z = expt_truncated_bidir_normal_from_qt( qt[i-1], qt[i], 0, 1, lg, lower);
      vz += z*z*den;
    }
  }
  double qu0 = qt[n-1];
  double ql0 = qt[0];
  if (lg) {
    qu0 = exp(qu0);
    ql0 = exp(ql0);
  }
  if (qu0 < 1.) {
    double tmp = -R::qnorm5( (1-qu0)*2, 0, 1, true, false);
    if ( tmp > 0 ) vz += (0.25 * R::pchisq(tmp*tmp, 3, false, false));
    else vz += (0.5 - 0.25 * R::pchisq(tmp*tmp, 3, false, false));
  }
  if (ql0 > 0.) {
    double tmp = -R::qnorm5( ql0*2, 0, 1, true, false);
    if ( tmp > 0 ) vz += (0.25 * R::pchisq(tmp*tmp, 3, false, false));
    else vz += (0.5 - 0.25 * R::pchisq(tmp*tmp, 3, false, false));
  }
  return vz;
}


// [[Rcpp::export]]
List binom_multi_bidir_n(IntegerVector xs, NumericVector sizes, NumericVector ps,
                         bool pos_only = true, bool var_adj = true, double approx_under = 1e-4) {
  int n = xs.size(); // n is the length of array
  if ( sizes.size() != n )
    stop("sizes must have same length to xs");
  if ( ps.size() != n )
    stop("ps must have same length to xs");
  NumericVector zs(n);
  NumericVector vz(n);
  NumericVector fadj(n);
  std::fill(fadj.begin(), fadj.end(), 1.);
  std::fill(vz.begin(), vz.end(), 0.);
  if (pos_only) {
    for (int i = 0; i < n; i++) {
      fadj[i] = 1. - R::dbinom(0, sizes[i], ps[i], false);
    }
  }
  double ql, qu;
  for (int i = 0; i < n; i++) {
    if (pos_only && xs[i] == 0) {
      zs[i] = 0; vz[i] = 0;
    } else {
      ql = R::pbinom(xs[i],sizes[i],ps[i],0,1) - log(fadj[i]);
      qu = (xs[i] == 0) ? ql : R::pbinom(xs[i]-1,sizes[i],ps[i],0,1) - log(fadj[i]);
      zs[i] = expt_truncated_bidir_normal_from_qt(ql, qu, 0, 1, 1);
    }
  }
  if (var_adj) {
    int k;
    double z,den, qu0, ql0, tmp;
    for (int i = 0; i < n; i++) {
      int median = R::qbinom(fadj[i]/2.0,sizes[i],ps[i],0,0);
      k = median;
      ql0 = R::pbinom(k,sizes[i],ps[i],0,1) - log(fadj[i]);
      qu0 = (k == 0) ? ql : R::pbinom(k-1,sizes[i],ps[i],0,1) - log(fadj[i]);
      z = expt_truncated_bidir_normal_from_qt( ql0, qu0, 0, 1, 1, 1);
      vz[i] += z * z * R::dbinom(k,sizes[i],ps[i],false)/fadj[i];
      qu = ql0;
      for (k = median + 1; k <= n; ++k) {
        den = R::dbinom(k,sizes[i],ps[i],false)/fadj[i];
        if (den < approx_under) break;
        ql = R::pbinom(k,sizes[i],ps[i],0,1) - log(fadj[i]);
        z = expt_truncated_bidir_normal_from_qt( ql, qu, 0, 1, 1, 1);
        vz[i] += z * z * den;
        qu = ql;
      }
      if (k <= sizes[i]) {
        tmp = -R::qnorm5((R::pbinom(k-1,sizes[i],ps[i],false,false))*2/fadj[i], 0, 1, true, false);
        if ( tmp > 0 ) vz[i] += (0.25 * R::pchisq(tmp*tmp, 3, false, false));
        else vz[i] += (0.5 - 0.25 * R::pchisq(tmp*tmp, 3, false, false));
      }
      ql = qu0;
      for (k = median - 1; k >= pos_only; --k) {
        den = R::dbinom(k,sizes[i],ps[i],false)/fadj[i];
        if (den < approx_under) break;
        qu = R::pbinom(k,sizes[i],ps[i],0,1) - log(fadj[i]);
        z = expt_truncated_bidir_normal_from_qt( ql, qu, 0, 1, 1, 1);
        vz[i] += z * z * den;
        ql = qu;
      }
      if (k >= pos_only) {
        tmp = -R::qnorm5((R::pbinom(k,sizes[i],ps[i],false,false))*2/fadj[i], 0, 1, true, false);
        if ( tmp > 0 ) vz[i] += (0.25 * R::pchisq(tmp*tmp, 3, false, false));
        else vz[i] += (0.5 - 0.25 * R::pchisq(tmp*tmp, 3, false, false));
      }
    }
  } else {
    std::fill(vz.begin(), vz.end(), 1.);
  }
  return List::create(Named("zs") = zs,
                      Named("variance") = vz);
}
