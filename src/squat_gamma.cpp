#include <Rcpp.h>
//#include <cmath>
using namespace Rcpp;

#define LOGDIFF(a,b) ((a > b) ? a + log(1.0-exp(b-a)) : b + log(1.0-exp(a-b))) // |a-b| in logscale
#define LOGADD(a,b) ((a < b) ? b + log(1.0+exp(a-b)) : a + log(1.0+exp(b-a))   // a+b in logscale
#define LOGTWO 0.6931471805599452862268 // log(2.0)
#define LOGSMALL -13.81551  // log(1e-6)
#define SMALL 1e-6;
#define MIN_SD 0.01 // minimum SD

//' @title expt_truncated_gamma_from_qt
//' @description
//' Calculate expectation from truncated Gamma
//'
//' @param a lower quantile
//' @param b upper quantile
//' @param alpha shape parameter of Gamma
//' @param beta rate parameter of Gamma
//' @param lg compute log-scale
//' @return A single expected z-score for directional test
// [[Rcpp::export]]
double expt_truncated_gamma_from_qt(double a, double b, double alpha, double beta, bool lg=0, bool lower=1) {
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
    if ( LOGDIFF(a,b) < LOGSMALL ) {return R::qgamma(a,alpha,1/beta,lower,lg);}
    denom = ( (b>a) ? exp(LOGDIFF(a,b)) : -exp(LOGDIFF(a,b)) );
  } else {
    if ( abs(a-b) < 1e-6 ) {
      return R::qgamma(a,alpha,1/beta,lower,lg);
    }
    denom = b-a;
  }
  return alpha/beta *
    ( R::pgamma(R::qgamma(b,alpha,1/beta,1,lg),alpha+1,1/beta,1,0) -
    R::pgamma(R::qgamma(a,alpha,1/beta,1,lg),alpha+1,1/beta,1,0) ) / denom;
}


//' @title expt_truncated_bidir_gamma_from_qt
//' @description
//' Calculate expectation from truncated Gamma - folded
//'
//' @param ql lower quantile
//' @param qu upper quantile
//' @param alpha shape parameter of Gamma
//' @param beta rate parameter of Gamma
//' @param lg compute log-scale
//' @return A single expected z-score for overdispersion test
// [[Rcpp::export]]
double expt_truncated_bidir_gamma_from_qt( double ql, double qu, double alpha, double beta, bool lg=0, bool lower=1 ) {
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
      ez = expt_truncated_gamma_from_qt( LOGDIFF(ql2,0), LOGDIFF(qu2,0), alpha, beta, lg, lower );
    } else if (qu <= log(0.5)) {
      ez = expt_truncated_gamma_from_qt( LOGDIFF(qu2,0), LOGDIFF(ql2,0), alpha, beta, lg, lower );
    } else {
      double w = ( 0.5-exp(ql) ) / ( exp(qu)-exp(ql) );
      ez = w * expt_truncated_gamma_from_qt( R_NegInf, LOGDIFF(ql2,0), alpha, beta, lg, lower ) +
        (1.-w) * expt_truncated_gamma_from_qt( R_NegInf, LOGDIFF(qu2,0), alpha, beta, lg, lower );
    }
  } else {
    if (!lower) {
      double tmp = 1-qu;
      qu = 1-ql;
      ql = tmp;
      lower=1;
    }
    if (ql > 0.5) {
      ez = expt_truncated_gamma_from_qt( 2.*ql-1., 2.*qu-1., alpha, beta, lg, lower );
    } else if (qu <= 0.5) {
      ez = expt_truncated_gamma_from_qt( 1.-2.*qu, 1.-2.*ql, alpha, beta, lg, lower );
    } else {
      double w = ( 0.5-ql ) / (qu - ql);
      ez = w * expt_truncated_gamma_from_qt(0, 1.-2.*ql, alpha, beta, lg, lower ) +
        (1.-w) * expt_truncated_gamma_from_qt(0, 2.*qu-1., alpha, beta, lg, lower );
    }
  }
  return ez;

}


//' @title squat_single_binom_unidir_g
//' @description
//' Function to calculate unidirectional z-score for single binomial distribution using Gamma approximation
//'
//' @param x            observed count
//' @param n            total number of trials
//' @param p            0-1 ranged binomial probability
//' @param alpha        shape parameter of the Gamma distribution
//' @param beta         rate parameter of the Gamma distribution
//' @param var_adj      perform variance adjustment if TRUE
//' @param approx_under threshold of binomial density to perform approximation during variance adjustment
//' @return A z-score and its standard deviation
// [[Rcpp::export]]
NumericVector squat_single_binom_unidir_g(int x, int n, double p, double alpha, double beta,
                                          bool var_adj = true,
                                          double approx_under = 1e-4, bool lower = true) {
  double sd = 1.0;  // without variance adjustment, assume unit variance
  double mu = n*p;
  if ( var_adj ) {  // when variance adjustment is enabled, calculate new variance
    int i;
    double sum = 0.0;
    double ez, f0, gf0, gf1, den;
    f0 = gf0 = 0.0;
    if ( approx_under == 0 ) { // calculate exact variance in log-scale to avoid underflow (could be slow)
      for(i=0; i <= n; ++i) {
        den = R::dbinom(i, n, p, false); // f(x)=Pr(x=i)
        f0 = f0 + den;
        gf1 = R::pgamma(R::qgamma(f0,alpha,1/beta,true,false),alpha+1,1/beta,1,0);
        ez = alpha / beta * (gf1 - gf0) / den;
        gf0 = gf1;
        if (!std::isnan(ez)) sum += den * pow( (ez-mu), 2 ); // otherwise, add E[z|x]^2f(x)
      }
    }
    else { // approximate variance with threshold
      double fnp, gfnp;
      sum = 0;
      int np = (int)(n*p); // median is somewhere between np and np+1
      fnp = f0 = R::pbinom(np, n, p, true, false);
      gfnp= gf0 = R::pgamma(R::qgamma(f0,alpha,1/beta,true,false),alpha+1,1/beta,1,0);
      for(i=np+1; i <= n; ++i) { // add upper-tail exact propobabilities
        den = R::dbinom(i, n, p, false); // f(x)=Pr(x=i)
        f0 = f0 + den;
        if (f0 > 1) {f0 = 1-1e-8;}
        gf1 = R::pgamma(R::qgamma(f0,alpha,1/beta,true,false),alpha+1,1/beta,1,0);
        ez = alpha / beta * (gf1 - gf0) / den;
        gf0 = gf1;
        if ( ( den < approx_under ) && ( i - 1 > n * p ) )
          break; // stop exact calculation and perform approximation
        sum += den * pow( (ez-mu), 2 ); // otherwise, add E[z|x]^2f(x)
      }
      if ( i <= n ) { // approximation needed
        sum += (alpha*(alpha+1)/beta/beta) * R::pgamma(f0,alpha+2,1/beta,false,false) -
          (2*alpha*alpha/beta/beta) * R::pgamma(f0,alpha+1,1/beta,false,false) +
          (alpha*alpha/beta/beta) * R::pgamma(f0,alpha,1/beta,false,false);
      }
      f0 = fnp; // f1=F(x=np)
      gf1= gfnp;
      for(i=np; i >= 0; --i) { // add lower-tail exact probabilities
        den = R::dbinom(i, n, p, false); // f(x)=Pr(x=i)
        f0 = f0 - den; // f0=F(x-1)=F(x)-f(x)
        if (f0 < 1e-8) {f0 = 0.0;}
        gf0= R::pgamma(R::qgamma(f0,alpha,1/beta,true,false),alpha+1,1/beta,1,0);
        ez = alpha / beta * (gf1 - gf0) / den;
        gf1 = gf0;
        if ( den < approx_under )
          break; // stop exact calculation and perform approximation
        sum += den * pow( (ez-mu), 2 ); // otherwise, add E[z|x]^2f(x)
      }
      if ( i >= 0 ) { // approximation needed
        sum += (alpha*(alpha+1)/beta/beta) * R::pgamma(f0,alpha+2,1/beta,true,false) -
          (2*alpha*alpha/beta/beta) * R::pgamma(f0,alpha+1,1/beta,true,false) +
          (alpha*alpha/beta/beta) * R::pgamma(f0,alpha,1/beta,true,false);
      }
    }
    sd = std::sqrt(sum); // update the adjusted variance
  }
  // calculate E[z|x]
  double logp1 = R::pbinom(x,n,p,0,1);
  double logp0 = R::pbinom(x-1,n,p,0,1);
  double Ez = expt_truncated_gamma_from_qt(logp1, logp0, alpha, beta, 1, lower);
  return NumericVector::create(Ez, sd);
}


//' @title squat_multi_binom_dir_g
//' @description
//' A function to generate two sample directional z scores based on exact quantiles from binomial distribution
//'
//' @param xs A integer vector containing the list of observed counts.
//' @param sizes A integer vector containg the list of total counts. Must be the same length with xs or a constant
//' @param ps A numeric vector containing the binomial probability for each observations. Must be the same length with ps or a constant
//' @param ys A binary vector indicating two groups.
//' @param var_adj Apply variance adjustment to improve power
//' @param approx_under Perform approximation in variance adjustment for Pr(X=x) smaller than the value
//' @return A list including z-scores and its two moments
// [[Rcpp::export]]
List squat_multi_binom_dir_g(IntegerVector xs, NumericVector sizes, NumericVector ps,
                             NumericVector ys, bool var_adj = true, double approx_under = 1e-4) {
  int n = xs.size(); // n is the length of array
  if ( ( sizes.size() != n ) && ( sizes.size() != 1 ) )
    stop("sizes must have same length to xs, or length of 1");
  if ( ( ps.size() != n ) && ( ps.size() != 1 ) )
    stop("ps must have same length to xs, or length of 1");
  if ( ys.size() != n )
    stop("ys must have same length to xs");
  int group_count = 0;
  for (int i=1; i<n; ++i) {
    if (ys[i] == 1) {
      group_count++;
    } else {
      ys[i] = 0;
    }
  }

  double beta = 0.0, alpha_sum = 0.0;
  NumericVector alpha(n);
  if (ps.size() == 1) {
    beta = 1./(1.-ps[0]);
    if (sizes.size() == 1) {
      for (int i=0; i<n; i++) {
        alpha[i] = sizes[0]*ps[0]*beta;
      }
      alpha_sum = n * alpha[0];
    } else {
      for (int i=0; i<n; i++) {
        alpha[i] = sizes[i]*ps[0]*beta;
        alpha_sum += alpha[i];
      }
    }
  } else {
    for (int i=0; i<n; i++) {
      beta += 1./(1.-ps[i]);
    }
    beta /= n;
    if (sizes.size() == 1) {
      for (int i=0; i<n; i++) {
        alpha[i] = sizes[0]*ps[i]*beta;
        alpha_sum += alpha[i];
      }
    } else {
      for (int i=0; i<n; i++) {
        alpha[i] = sizes[i]*ps[i]*beta;
        alpha_sum += alpha[i];
      }
    }
  }

  NumericVector zs(n);
  NumericVector res(2);
  double v_z = 0.0;
  for(int i=0; i < n; ++i) {
    res = squat_single_binom_unidir_g(xs[i], (sizes.size() == 1) ? sizes[0] : sizes[i],
                                      (ps.size() == 1) ? ps[0] : ps[i],
                                      alpha[i], beta, var_adj, approx_under, (1-ys[i]));
    zs[i] = res[0];
    v_z += res[1]*res[1];
  }
  return List::create(Named("zs") = zs,
                      Named("mean") = alpha_sum / beta,
                      Named("variance") =  v_z);
}
