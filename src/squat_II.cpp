#include <Rcpp.h>
#include "shared.h"
// //#include <cmath>
using namespace Rcpp;

//' Expectation of truncated normal from quantiles
//' Calculate expectation from truncated Normal
//' 
//' @param a lower quantile
//' @param b upper quantile
//' @param mu mean of the Normal
//' @param sd standard deviation of the Normal
//' @param lg compute log-scale
//' @param lower use lower tial quantiles
//' @return A single expected z-score for directional test
// [[Rcpp::export]]
double expt_truncated_normal_from_qt(double a, double b, 
                                     double mu=0, double sd=1, 
                                     bool lg=0, bool lower=1) {
  double denom;
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
    if ( std::abs(b-a) < 1e-6 ) {
      return R::qnorm(b,mu,sd,lower,lg);
    }
    denom = b-a;
  }
  return mu + sd * (R::dnorm(R::qnorm(a,mu,sd,1,lg),mu,sd,0) -
                    R::dnorm(R::qnorm(b,mu,sd,1,lg),mu,sd,0) ) / denom;
}


//' Expectation of truncated normal from bidirectional quantiles
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


//' Compute multiple bi-directional expected z-scores from quantiles 
//' @description
//' A function to generate bi-directional z scores based on input quantiles (from an arbitrary distribution)
//'
//' @param ql A vector of lower quantiles
//' @param qu A vector of upper quantiles
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


//' Variance of bi-directional expected z-scores from quantiles
//' @name bidir_etn_var_from_qt
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
    if (lg && (qt[n-1] < log(0.5) || qt[0] > log(0.5)) ) {
      stop("The exact calculation interval must cover 0.5.");
    }
    if (!lg && (qt[n-1] < 0.5 || qt[0] > 0.5) ) {
      stop("The exact calculation interval must cover 0.5.");
    }
    double z;
    for (int i = 1; i < n; i++) {
      double den = lg ? exp(LOGDIFF(qt[i-1],qt[i])) : (qt[i]-qt[i-1]);
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


//' Binomial bi-directional test 
//' @name binom_multi_bidir_n
//' 
//' @param xs A vector of non-negative counts representing observed data
//' @param sizes A vector of positive values representing total counts
//' @param ps A non-negative vector of probabilities
//' @param pos_only Ignore zero observations
//' @param var_adj Adjust variance to improved power
//' @param approx_under Use approximate variance calculation when Pr(X)<value
//' @return A list including z-scores and their variance
// [[Rcpp::export]]
List binom_multi_bidir_n(IntegerVector xs, IntegerVector sizes, NumericVector ps,
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



//' Beta-binomial bi-directional test 
//' @name betabinom_multi_bidir_n
//' 
//' @param xs A vector of non-negative counts representing observed data
//' @param sizes A vector of positive values representing total counts
//' @param alphas A non-negative vector of parameters of the beta distribution
//' @param betas A non-negative vector of parameters of the beta distribution
//' @param pos_only Ignore zero observations
//' @param var_adj Adjust variance to improved power
//' @param approx_under Use approximate variance calculation when Pr(X)<value
//' @return A list including z-scores and their variance
// [[Rcpp::export]]
List betabinom_multi_bidir_n(NumericVector xs, NumericVector sizes,
                             NumericVector alphas, NumericVector betas,
                             bool pos_only = true, bool var_adj = true,
                             double approx_under = 1e-4) {
  // double log_thres = log(approx_under);
  int n = xs.size(); // n is the length of array
  if ( sizes.size() != n )
    stop("sizes must have same length to xs");
  if ( alphas.size() != n )
    stop("alphas must have same length to xs");
  NumericVector zs(n);
  NumericVector vz(n);
  NumericVector fadj(n);
  NumericVector fadjorg(n);
  std::fill(fadj.begin(), fadj.end(), 0.);
  std::fill(fadjorg.begin(), fadjorg.end(), 1.);
  std::fill(vz.begin(), vz.end(), 0.);
  if (pos_only) {
    fadj = cpp_pbbinom(fadj, sizes, alphas, betas, false, true);
    fadjorg = exp(fadj);
  }
  NumericVector qls = cpp_pbbinom(xs,sizes,alphas,betas,0,0);
  NumericVector qus = qls + cpp_dbbinom(xs,sizes,alphas,betas,0);
  qls = log(qls / fadjorg);
  qus = log(qus / fadjorg);
  zs = multi_bidir_zs_from_qt_n(qls, qus, 1, 1);  
  
  double ql, qu;
  if (var_adj) {
    int k, minx, maxx;
    double z, den, qu0, ql0, tmp;
    for (int i = 0; i < n; i++) {
      int median = sizes[i] * alphas[i] / (alphas[i] + betas[i]);
      if (median <= 0) {median=1;}
      k = median;
      den = cpp_dbbinom_one(k,sizes[i],alphas[i],betas[i],false)/fadjorg[i];
      ql0 = cpp_pbbinom_one(k,sizes[i],alphas[i],betas[i],0,1) - fadj[i];
      qu0 = (k == 0) ? ql0 : cpp_pbbinom_one(k-1,sizes[i],alphas[i],betas[i],0,1) - fadj[i];
      z = expt_truncated_bidir_normal_from_qt( ql0, qu0, 0, 1, 1, 1);
      vz[i] += z * z * den;
      qu = ql0;
      for (k = median + 1; k <= n; ++k) {
        den = cpp_dbbinom_one(k,sizes[i],alphas[i],betas[i],false)/fadjorg[i];
        if (den < approx_under) break;
        ql = cpp_pbbinom_one(k,sizes[i],alphas[i],betas[i],0,1) - fadj[i];
        z = expt_truncated_bidir_normal_from_qt( ql, qu, 0, 1, 1, 1);
        vz[i] += z * z * den;
        qu = ql;
      }
      maxx = k;
      if (k <= sizes[i]) {
        tmp = -R::qnorm5((cpp_pbbinom_one(k-1,sizes[i],alphas[i],betas[i],false,false))*2/fadjorg[i], 0, 1, true, false);
        if ( tmp > 0 ) vz[i] += (0.25 * R::pchisq(tmp*tmp, 3, false, false));
        else vz[i] += (0.5 - 0.25 * R::pchisq(tmp*tmp, 3, false, false));
      }
      ql = qu0;
      for (k = median - 1; k >= pos_only; --k) {
        den = cpp_dbbinom_one(k,sizes[i],alphas[i],betas[i],false)/fadjorg[i];
        if (den < approx_under) break;
          qu = cpp_pbbinom_one(k,sizes[i],alphas[i],betas[i],0,1) - fadj[i];
          z = expt_truncated_bidir_normal_from_qt( ql, qu, 0, 1, 1, 1);
          vz[i] += z * z * den;
          ql = qu;
      }
      if (k < 0) minx = 0; 
      else  minx = k;
      if (k >= pos_only) {
        tmp = -R::qnorm5((cpp_pbbinom_one(k,sizes[i],alphas[i],betas[i],true,false))*2/fadjorg[i], 0, 1, true, false);
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