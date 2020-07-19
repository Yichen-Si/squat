#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#define LOGDIFF(a,b) ((a > b) ? a + log(1.0-exp(b-a)) : b + log(1.0-exp(a-b))) // |a-b| in logscale
#define LOGADD(a,b) ((a < b) ? b + log(1.0+exp(a-b)) : a + log(1.0+exp(b-a))   // a+b in logscale
#define LOG2 0.6931471805599452862268 // log(2.0)
#define LOGSMALL -13.81551  // log(1e-6)
#define MIN_SD 0.01 // minimum SD


//' Internal function to calculate \phi(\Phi^{-1}(F(x)))
//' @param x observed count
//' @param n total number of trials
//' @param p 0-1 range binomial probability
//' @param lg compute log-scale
//' @return \phi(\Phi^{-1}(F(x))) where F(x) = pbinom(x,n,p)
double dnorm_qnorm_pbinom(int x, int n, double p, bool lg) {
  if ( ( x < 0 ) || ( x == n ) ) 
    return lg ? R_NegInf : 0;
  else {
    bool lt = x < n*p; // compute lower tail
    return R::dnorm4(R::qnorm5(R::pbinom(x,n,p,lt,lg),0,1,lt,lg),0,1,lg);
  }
}

//' Internal function to calculate \phi(\Phi^{-1}(2*F(x))) or \phi(\Phi^{-1}(2(F(x)-F(0))/(1-F(0))))
//' @param x observed count (smaller than median)
//' @param n total number of trials
//' @param p 0-1 range binomial probability
//' @param lg compute log-scale
//' @param pos_only consider positive values only
//' @param f0 Pr(X=0) or log Pr(X=0)
//' @return \phi(\Phi^{-1}(2*F(x))) (pos_only=false) or \phi(\Phi^{-1}(2(F(x)-F(0))/(1-F(0)))) (pos_only=true)
double dnorm_qnorm_2pbinom_lt(int x, int n, double p, bool lg, bool pos_only, double f0) {
  //Rprintf("x=%d, n=%d, p=%lg lg=%d pos_only=%d, f0=%lg\n",x,n,p,lg,pos_only,f0);
  if ( x < pos_only )
    return lg ? R_NegInf : 0;
  else if ( lg ) { // log scale
    double logf = R::pbinom(x,n,p,true,true); // log Pr(X<=x)
    if ( pos_only ) { // use f0
      return R::dnorm4(R::qnorm5(logf+log(1.0-exp(f0-logf))+LOG2-log(1-exp(f0)),0,1,true,true),0,1,true);
    }
    else 
      return R::dnorm4(R::qnorm5(logf+LOG2,0,1,true,true),0,1,true);
  }
  else { // linear scale
    if ( pos_only )
      return R::dnorm4(R::qnorm5((R::pbinom(x,n,p,true,false)-f0)*2/(1-f0),0,1,true,false),0,1,false);
    else
      return R::dnorm4(R::qnorm5(R::pbinom(x,n,p,true,false)*2,0,1,true,false),0,1,false);
  }
}

//' Internal function to calculate \phi(\Phi^{-1}(2*F_c(x))) or \phi(\Phi^{-1}(2(F_c(x))/(1-F(0))))
//' @param x observed count (greater than median)
//' @param n total number of trials
//' @param p 0-1 range binomial probability
//' @param lg compute log-scale
//' @param pos_only consider positive values only
//' @param f0 Pr(X=0) or log Pr(X=0)
//' @return \phi(\Phi^{-1}(2*F_c(x))) (pos_only=false) or \phi(\Phi^{-1}(2(F_c(x))/(1-F(0)))) (pos_only=true)
double dnorm_qnorm_2pbinom_ut(int x, int n, double p, bool lg, bool pos_only, double f0) {
  //Rprintf("x=%d, n=%d, p=%lg lg=%d pos_only=%d, f0=%lg\n",x,n,p,lg,pos_only,f0);
  if ( ( x == n ) || ( x < pos_only ) )
    return lg ? R_NegInf : 0;
  else if ( lg ) { // log scale
    if ( pos_only ) { // use f0
      return R::dnorm4(R::qnorm5(R::pbinom(x,n,p,false,true)+LOG2-log(1-exp(f0)),0,1,false,true),0,1,true);
    }
    else 
      return R::dnorm4(R::qnorm5(R::pbinom(x,n,p,false,true)+LOG2,0,1,false,true),0,1,true);
  }
  else { // linear scale
    if ( pos_only )
      return R::dnorm4(R::qnorm5(R::pbinom(x,n,p,false,false)*2/(1-f0),0,1,false,false),0,1,false);
    else
      return R::dnorm4(R::qnorm5(R::pbinom(x,n,p,false,false)*2,0,1,false,false),0,1,false);
  }
}

//' Internal function to calculate \phi(\Phi^{-1}(q)
//' @param q quantile value
//' @param lt quantile is lower tail
//' @param lg log-scale computation
//' @return \phi(\Phi^{-1}(q)
double dnorm_qnorm(double q, bool lt, bool lg) {
  return R::dnorm4(R::qnorm5(q,0,1,lt,lg),0,1,lg);
}

//' Function to calculate unidirectional z-score for single binomial distribution
//' 
//' @param x            observed count
//' @param n            total number of trials
//' @param p            0-1 ranged binomial probability
//' @param var_adj      perform variance adjustment if TRUE
//' @param approx_under threshold of binomial density to perform approximation during variance adjustment
// [[Rcpp::export]]
double squat_single_binom_unidir(int x, int n, double p, bool var_adj = true, double approx_under = 1e-4) {
  double sd = 1.0;  // without variance adjustment, assume unit variance
  if ( var_adj ) {  // when variance adjustment is enabled, calculate new variance
    int i;
    double sum = 0.0;
    if ( approx_under == 0 ) { // calculate exact variance in log-scale to avoid underflow (could be slow)
      double logden, logd1, logd0, logdiff, toadd;
      for(i=0; i <= n; ++i) {
        logden = R::dbinom(i, n, p, true);
        logd0 = (i > 0) ? logd1 : dnorm_qnorm_pbinom(i-1,n,p,true);
        logd1 = dnorm_qnorm_pbinom(i,n,p,true);
        logdiff = LOGDIFF(logd1,logd0);
        toadd = exp(2*logdiff - logden);   // sum is the variance
        if ( !isnan(toadd) ) sum += toadd; // to avoid the pbinom bug causing underflow occasionally
      }
    }
    else { // approximate variance in linear scale with threshold
      double fnp, dnp, f0, f1, d0, d1, den, tmp;
      sum = 0;
      int np = (int)(n*p); // median is somewhere between np and np+1
      fnp = f0 = R::pbinom(np, n, p, true, false); // F(x=np) 
      dnp = d0 = dnorm_qnorm(f0, false, false);    // \phi(\Phi^{-1}(F(np)))
      for(i=np+1; i <= n; ++i) { // add upper-tail exact propobabilities
        den = R::dbinom(i, n, p, false); // f(x)=Pr(x=i)
        f1 = f0 + den; // f1=F(x)
        d1 = (i == n) ? 0 : dnorm_qnorm(f1,false,false); // d1=\phi(\Phi^{-1}(F(x)))
        if ( ( den < approx_under ) && ( i - 1 > n * p ) ) 
          break; // stop exact calculation and perform approximation
        sum += ( ( d1 - d0 ) * ( d1 - d0 ) / den ); // otherwise, add E[z|x]^2f(x)
        d0 = d1;
        f0 = f1;
      }
      if ( i <= n ) { // approximation needed
        tmp = -R::qnorm5(f0, 0, 1, false, false); // calculate z corresponding F(x-1)
        sum += (0.5 * R::pchisq(tmp*tmp, 3, false, false)); // integrate upper-tail assuming continuity
      }
      f1 = fnp; // f1=F(x=np)
      d1 = dnp; // d1=\phi(\Phi^{-1}(f1))
      for(i=np; i >= 0; --i) { // add lower-tail exact probabilities
        den = R::dbinom(i, n, p, false); // f(x)=Pr(x=i)
        f0 = f1 - den; // f0=F(x-1)=F(x)-f(x)
        d0 = (i == 0) ? 0 : dnorm_qnorm(f0,true,false); // d0=\phi(\Phi^{-1}(f1))
        if ( den < approx_under )
          break; // stop exact calculation and perform approximation
        sum += ( ( d1 - d0 ) * ( d1 - d0 ) / den ); // otherwise, add E[z|x]^2f(x)
        d1 = d0;
        f1 = f0;
      }
      if ( i >= 0 ) { // approximation needed
        tmp = R::qnorm5(f1, 0, 1, true, false); // calculate z corresponding F(x-1)
        sum += (0.5 * R::pchisq(tmp*tmp, 3, false, false)); // integrate lower-tail assuming continuity
      }
    }
    sd = std::sqrt(sum); // update the adjusted variance
    //Rprintf("sd = %lg\n", sd);
  }
  // calculate E[z|x]
  double Ez, d1, d0, logden;
  logden = R::dbinom(x,n,p,true);
  if ( logden < LOGSMALL ) { // density is less than 1e-6 - calculate in log-scale for precision
    d0 = dnorm_qnorm_pbinom(x-1,n,p,true); // d0 = \phi(\Phi^{-1}(F(x-1)))
    d1 = dnorm_qnorm_pbinom(x,n,p,true);   // d1 = \phi(\Phi^{-1}(F(x)))
    // Ez = (d0-d1)/den calculated in log-scale
    Ez = (d0 < d1) ? -exp( d1 + log(1.0-exp(d0-d1)) - logden ) : exp( d0 + log(1.0-exp(d1-d0)) - logden );
  }
  else {
    d0 = dnorm_qnorm_pbinom(x-1,n,p,false); // d0 = \phi(\Phi^{-1}(F(x-1)))
    d1 = dnorm_qnorm_pbinom(x,n,p,false);   // d1 = \phi(\Phi^{-1}(F(x)))
    Ez = (d0-d1)/exp(logden); // Ez = (d0-d1)/den calculated in linear scale
  }
  if ( sd < MIN_SD ) sd = MIN_SD; // prevent sd from becoming zero or too small
  Ez /= sd;
  return Ez;
}

//' Function to calculate bidirectional z-score for single binomial distribution
//' 
//' @param x            observed count
//' @param n            total number of trials
//' @param p            0-1 ranged binomial probability
//' @param pos_only     ignore zeros if TRUE (x must be positive)
//' @param var_adj      perform variance adjustment if TRUE
//' @param approx_under threshold of binomial density to perform approximation during variance adjustment
// [[Rcpp::export]]
double squat_single_binom_bidir(int x, int n, double p, bool pos_only = true, bool var_adj = true, double approx_under = 1e-4) {
  double sd = 1.0;  // without variance adjustment assume unit variance
  double logden0, den0, pdenom, logpdenom, logden, den, d0, d1, Ez, tmp, sum;
  int median, i;
  
  logden0 = pos_only ? R::dbinom(0,n,p,true) : R_NegInf; // log Pr(X=0) when pos_only
  den0 = pos_only ? exp(logden0) : 0; // Pr(X=0) when pos_only
  pdenom = pos_only ? 1-den0 : 1.; // Pr(X>0) when pos_only 
  logpdenom = log(pdenom); // log Pr(X>0) when pos_only
  median = R::qbinom(pdenom/2.0,n,p,false,false); // median among the valid range
  
  if ( var_adj ) {  // perform adjustment of variance
    if ( approx_under == 0 ) { // exact calculation of variance
      logden = R::dbinom(median, n, p, true); // log density
      d0 = dnorm_qnorm_2pbinom_lt(median-1,n,p,true,pos_only,logden0); 
      d1 = dnorm_qnorm_2pbinom_ut(median,n,p,true,pos_only,logden0);
      Ez = (d0 < d1) ? -exp( d1 + log(1.0+exp(d0-d1)) - logden - LOG2 + logpdenom ) : -exp( d0 + log(1.0+exp(d1-d0)) - logden - LOG2 + logpdenom);
      sum = (Ez*Ez*exp(logden-logpdenom));
      for(i=median+1; i <= n; ++i) {
        logden = R::dbinom(i, n, p, true); // log density
        d0 = dnorm_qnorm_2pbinom_ut(i-1,n,p,true,pos_only,logden0);
        d1 = dnorm_qnorm_2pbinom_ut(i,n,p,true,pos_only,logden0);
        Ez = (d0 < d1) ? -exp( d1 + log(1.0-exp(d0-d1)) - logden -LOG2 + logpdenom ) : exp( d0 + log(1.0-exp(d1-d0)) - logden - LOG2 + logpdenom);
        if ( !isnan(Ez) ) sum += (Ez*Ez*exp(logden-logpdenom)); // avoid occasional pbinom underflow 
      }
      for(i=median-1; i >= pos_only; --i) {
        logden = R::dbinom(i, n, p, true); // log density
        d0 = dnorm_qnorm_2pbinom_lt(i-1,n,p,true,pos_only,logden0);
        d1 = dnorm_qnorm_2pbinom_lt(i,n,p,true,pos_only,logden0);
        Ez = (d0 < d1) ? exp( d1 + log(1.0-exp(d0-d1)) - logden - LOG2 + logpdenom ) : -exp( d0 + log(1.0-exp(d1-d0)) - logden - LOG2 + logpdenom);
        if ( !isnan(Ez) ) sum += (Ez*Ez*exp(logden-logpdenom)); // avoid occasional pbinom underflow
      }
    }
    else { // perform approximation when calculating variance
      den = R::dbinom(median, n, p, false); // log density
      d0 = dnorm_qnorm_2pbinom_lt(median-1,n,p,false,pos_only,den0);
      d1 = dnorm_qnorm_2pbinom_ut(median,n,p,false,pos_only,den0);
      sum = ((d0+d1)*(d0+d1)/den/4*pdenom); // always calculate exact probability for median
      for(i=median+1; i <= n; ++i) { // upper tail probability
        den = R::dbinom(i, n, p, false);
        if ( den < approx_under ) break; // start approximation if density is small
        d0 = dnorm_qnorm_2pbinom_ut(i-1,n,p,false,pos_only,den0);
        d1 = dnorm_qnorm_2pbinom_ut(i,n,p,false,pos_only,den0);
        sum += ((d0-d1)*(d0-d1)/den/4*pdenom);
      }
      if ( i <= n ) { // approximation is needed
        tmp = -R::qnorm5((R::pbinom(i-1,n,p,false,false))*2/pdenom, 0, 1, true, false);
        if ( tmp > 0 ) sum += (0.25 * R::pchisq(tmp*tmp, 3, false, false));
        else sum += (0.5 - 0.25 * R::pchisq(tmp*tmp, 3, false, false)); 
      }
      for(i=median-1; i >= pos_only; --i) {
        den = R::dbinom(i, n, p, false);
        if ( den < approx_under ) break;
        d0 = dnorm_qnorm_2pbinom_lt(i-1,n,p,false,pos_only,den0);
        d1 = dnorm_qnorm_2pbinom_lt(i,n,p,false,pos_only,den0);
        sum += ((d0-d1)*(d0-d1)/den/4*pdenom);
      }
      if ( i >= pos_only ) {
        tmp = -R::qnorm5((R::pbinom(i,n,p,true,false)-den0)*2/pdenom, 0, 1, true, false);
        if ( tmp > 0 ) sum += (0.25 * R::pchisq(tmp*tmp, 3, false, false));
        else sum += (0.5 - 0.25 * R::pchisq(tmp*tmp, 3, false, false));
      }
    }
    sd = std::sqrt(sum);
    //Rprintf("sd=%lg\n",sd);
  }
  
  // calculate E[z|x]
  logden = R::dbinom(x, n, p, true); // log density
  if ( logden < LOGSMALL ) { // density < 1e-6, compute in log-scale
    if ( x == median ) {
      d0 = dnorm_qnorm_2pbinom_lt(x-1,n,p,true,pos_only,logden0);
      d1 = dnorm_qnorm_2pbinom_ut(x,n,p,true,pos_only,logden0);
      Ez = (d0 < d1) ? 
            -exp( d1 + log(1.0+exp(d0-d1)) - logden - LOG2 + logpdenom) : 
            -exp( d0 + log(1.0+exp(d1-d0)) - logden - LOG2+ logpdenom);
    }    
    else if ( x < median ) {
      d0 = dnorm_qnorm_2pbinom_lt(x-1,n,p,true,pos_only,logden0);
      d1 = dnorm_qnorm_2pbinom_lt(x,n,p,true,pos_only,logden0);
      Ez = (d0 < d1) ? 
            exp( d1 + log(1.0-exp(d0-d1)) - logden - LOG2 + logpdenom ) : 
            -exp( d0 + log(1.0-exp(d1-d0)) - logden - LOG2 + logpdenom);
    }
    else {
      d0 = dnorm_qnorm_2pbinom_ut(x-1,n,p,true,pos_only,logden0);
      d1 = dnorm_qnorm_2pbinom_ut(x,n,p,true,pos_only,logden0);
      Ez = (d0 < d1) ? 
            -exp( d1 + log(1.0-exp(d0-d1)) - logden -LOG2 + logpdenom) : 
            exp( d0 + log(1.0-exp(d1-d0)) - logden - LOG2 + logpdenom);
    }
  }
  else {
    if ( x == median ) {
      d0 = dnorm_qnorm_2pbinom_lt(x-1,n,p,false,pos_only,den0);
      d1 = dnorm_qnorm_2pbinom_ut(x,n,p,false,pos_only,den0);
      Ez = -(d0+d1)/exp(logden)/2*pdenom;
    }
    else if ( x < median ) {
      d0 = dnorm_qnorm_2pbinom_lt(x-1,n,p,false,pos_only,den0);
      d1 = dnorm_qnorm_2pbinom_lt(x,n,p,false,pos_only,den0);
      Ez = (d1-d0)/exp(logden)/2*pdenom;
    }
    else { // x > median
      d0 = dnorm_qnorm_2pbinom_ut(x-1,n,p,false,pos_only,den0);
      d1 = dnorm_qnorm_2pbinom_ut(x,n,p,false,pos_only,den0);
      Ez = (d0-d1)/exp(logden)/2*pdenom;
    }
  }
  if ( sd < MIN_SD ) sd = MIN_SD;
  return Ez/sd;
}

//' A function to generate unidirectional z scores based on exact quantiles from binomial distribution
//' 
//' @param xs A integer vector containing the list of observed counts. 
//' @param sizes A integer vector containg the list of total counts. Must be the same length with xs or a constant
//' @param ps A numeric vector containing the binomial probability for each observations. Must be the same length with ps or a constant
//' @param var_adj Apply variance adjustment to improve power
//' @param approx_under Perform approximation in variance adjustment for Pr(X=x) smaller than the value
//' @return A vector of z-scores corresponding to expected overdispersion z-scores from squat
// [[Rcpp::export]]
NumericVector squat_multi_binom_unidir (IntegerVector xs, NumericVector sizes, NumericVector ps, bool var_adj = true, double approx_under = 1e-4) {
  int n = xs.size(); // n is the length of array
  if ( ( sizes.size() != n ) && ( sizes.size() != 1 ) )
    stop("sizes must have same length to xs, or length of 1");
  if ( ( ps.size() != n ) && ( ps.size() != 1 ) )
    stop("ps must have same length to xs, or length of 1");
  NumericVector zs(n);
  for(int i=0; i < n; ++i) {
    zs[i] = squat_single_binom_unidir(xs[i], (sizes.size() == 1) ? sizes[0] : sizes[i], 
                                      (ps.size() == 1) ? ps[0] : ps[i], var_adj, approx_under);
  }
  return zs;
}

//' A function to generate bidirectional z scores based on exact quantiles from binomial distribution
//' 
//' @param xs A integer vector containing the list of observed counts. 
//' @param sizes A integer vector containg the list of total counts. Must be the same length with xs or a constant
//' @param ps A numeric vector containing the binomial probability for each observations. Must be the same length with ps or a constant
//' @param pos_only Ignore zeros in the distribution
//' @param var_adj Apply variance adjustment to improve power
//' @param approx_under Perform approximation in variance adjustment for Pr(X=x) smaller than the value
//' @return A vector of z-scores corresponding to expected overdispersion z-scores from squat
// [[Rcpp::export]]
NumericVector squat_multi_binom_bidir (IntegerVector xs, NumericVector sizes, NumericVector ps, bool pos_only = true, bool var_adj = true, double approx_under = 1e-4) {
  int n = xs.size(); // n is the length of array
  if ( ( sizes.size() != n ) && ( sizes.size() != 1 ) )
    stop("sizes must have same length to xs, or length of 1");
  if ( ( ps.size() != n ) && ( ps.size() != 1 ) )
    stop("ps must have same length to xs, or length of 1");
  NumericVector zs(n);
  for(int i=0; i < n; ++i) {
    if ( pos_only && ( xs[i] == 0 ) ) {
      zs[i] = std::numeric_limits<double>::quiet_NaN(); 
    }
    else {
      zs[i] = squat_single_binom_bidir(xs[i], sizes.size() == 1 ? sizes[0] : sizes[i], 
                                      ps.size() == 1 ? ps[0] : ps[i], pos_only, var_adj, approx_under);
    }
  }
  return zs;
}
  

//' A function to generate directional z scores based on exact quantiles from binomial distribution
//' 
//' @param xs A integer vector containing the list of observed counts. Must be positive if posOnly = TRUE.
//' @param sizes A integer vector containg the list of total counts. Must be the same length with xs
//' @param ps A numeric vector containing the binomial probability for each observations. Must be the same length with ps
//' @return A vector of z-scores corresponding to expected overdispersion z-scores from squat
// [[Rcpp::export]]
NumericVector squat_z_bd (IntegerVector xs, NumericVector sizes, NumericVector ps) {
  int n = xs.size(); // n is the length of array
  if ( (n != sizes.size()) || ( n != ps.size() ) ) 
    stop("xs, sizes, ps must have the same size");
  NumericVector zs(n);
  double ld, lpl, lpu, ll, lu;
  for(int i=0; i < n; ++i) {
    if ( ( xs[i] < 0 ) || ( xs[i] > sizes[i] ) )
      stop("xs must be nonnegative and no greater than ns");
    if ( ( ps[i] <= 0 ) || ( ps[i] >= 1 ) )
      stop("ps must be a probability between 0 and 1");
    
    ld = R::dbinom(xs[i], sizes[i], ps[i], true); // log Pr(X=x)
    lpl = (xs[i] == 0)     ? R_NegInf : R::pbinom(xs[i]-1, sizes[i], ps[i], true, true); // log Pr(X<x)
    if ( std::isinf(lpl) ) lpl = R_NegInf;
    lpu = (xs[i] == sizes[i]) ? R_NegInf : R::pbinom(xs[i], sizes[i], ps[i], false, true);  // log [ Pr(X>x)/Pr(X>0) ]
    if ( std::isinf(lpu) ) lpu = R_NegInf;
    
    if ( lpl < lpu ) {  // lower quantile is smaller
      ll = R::dnorm( R::qnorm(lpl , 0, 1, true, true), 0, 1, false );
      lu = R::dnorm( R::qnorm(ld + log( 1 + exp ( lpl - ld) ), 0, 1, true, true), 0, 1, false );
      zs[i] = (ll - lu) / exp(ld);
    }
    else { // upper quantile is smaller
      ll = R::dnorm( R::qnorm(lpu , 0, 1, true, true), 0, 1, false );
      lu = R::dnorm( R::qnorm(ld + log( 1 + exp ( lpu - ld) ), 0, 1, true, true), 0, 1, false );
      zs[i] = (lu - ll) / exp(ld);
    }
  }
  return(zs);
} 

//' A function to generate overdispersion z scores based on exact quantiles from binomial distribution
//' 
//' @param xs A integer vector containing the list of observed counts. Must be positive if posOnly = TRUE.
//' @param sizes A integer vector containg the list of total counts. Must be the same length with xs
//' @param ps A numeric vector containing the binomial probability for each observations. Must be the same length with ps
//' @param posOnly Ignore zero observations and only consider positive values. All values of xs must be positive.
//' @return A vector of z-scores corresponding to expected overdispersion z-scores from squat
// [[Rcpp::export]]
NumericVector squat_z_bo (IntegerVector xs, NumericVector sizes, NumericVector ps, bool posOnly) {
  double log2 = 0.6931471805599452862268;
  double logThres = -46.0517; // log(1e-20);
  int n = xs.size();
  if ( (n != sizes.size()) || ( n != ps.size() ) ) 
    stop("xs, sizes, ps must have the same size");
  NumericVector zs(n);
  double l0, ld, lpl, lpu, ll, lu;
  for(int i=0; i < n; ++i) {
    if ( ( xs[i] < 0 ) || ( xs[i] > sizes[i] ) )
      stop("xs must be nonnegative and no greater than ns");
    if ( ( ps[i] <= 0 ) || ( ps[i] >= 1 ) )
      stop("ps must be a probability between 0 and 1");
    if ( posOnly ) { // positive only
      if ( xs[i] == 0 ) stop("Only positive values of xs are expected");
      l0 = log(1.0 - R::dbinom(0,sizes[i],ps[i],false));
      //l0 = R::pbinom(0, sizes[i], ps[i], false, true); // log Pr(X>0)
      //if ( std::isnan(l0) ) l0 = 0;
      ld = R::dbinom(xs[i], sizes[i], ps[i], true) - l0; // log [ Pr(X=x)/Pr(X>0) ]
      if ( xs[i] == 1 ) lpl =  R_NegInf;
      else {
        lpl = R::pbinom(xs[i]-1, sizes[i], ps[i], true, true); // log [ Pr(X<x)/Pr(X>0) ] 
        if ( std::isinf(lpl) ) lpl = R_NegInf;
        //Rprintf("%lf %lf %lf\n", l0, ld, lpl);
        else lpl += ( log(1 - exp(R::dbinom(0, sizes[i], ps[i], true) - lpl)) - l0 );
        //Rprintf("%lf %lf %lf\n", l0, ld, lpl);
        //lpl = R::pbinom(xs[i]-1, sizes[i], ps[i], true, true) - l0;
      }
      lpu = ((xs[i] == sizes[i]) ? R_NegInf : R::pbinom(xs[i], sizes[i], ps[i], false, true)) - l0;  // log [ Pr(X>x)/Pr(X>0) ]
      if ( std::isinf(lpu) ) lpu = R_NegInf;
    }
    else {
      ld = R::dbinom(xs[i], sizes[i], ps[i], true); // log Pr(X=x)
      lpl = (xs[i] == 0)     ? R_NegInf : R::pbinom(xs[i]-1, sizes[i], ps[i], true, true); // log Pr(X<x)
      if ( std::isinf(lpl) ) lpl = R_NegInf;
      lpu = (xs[i] == sizes[i]) ? R_NegInf : R::pbinom(xs[i],   sizes[i], ps[i], false, true);  // log [ Pr(X>x)/Pr(X>0) ]
      if ( std::isinf(lpu) ) lpu = R_NegInf;
    }
    
    if ( lpu > -log2 ) { // Q_l < Q_u < 0.5
      if ( ld < logThres ) {  // Use linear interpolation (avg between Q_l and Q_u) 
        // log(2Q_l + delta) = log(delta) + log(1 + exp(log2 + log Q_l - log delta))
        zs[i] = R::qnorm(ld + log(1.0 + exp(log2 + lpl - ld)), 0, 1, false, true);
      }
      else {
        //  (dnorm(qnorm(2*qu)) - dnorm(qnorm(2*ql)))/2/(qu-ql) 
        ll = R::dnorm( R::qnorm(lpl + log2, 0, 1, true, true), 0, 1, false );
        lu = R::dnorm( R::qnorm(ld + log( 1 + exp ( lpl - ld) ) + log2, 0, 1, true, true), 0, 1, false );
        zs[i] = (lu - ll) / 2 / exp(ld);
      }
    }
    else if ( lpl > -log2 ) { // Q_u > Q_l > 0.5
      if ( ld < logThres ) {  // Use linear interpolation (avg between Q_l and Q_u) 
        // log(2Q_l + delta) = log(delta) + log(1 + exp(log2 + log Q_l - log delta))
        zs[i] = R::qnorm(ld + log(1.0 + exp(log2 + lpu - ld)), 0, 1, false, true);
      }
      else {
        // (dnorm(qnorm(2*(1-ql))) - dnorm(qnorm(2*(1-qu))))/2/(qu-ql)
        ll = R::dnorm( R::qnorm(lpu + log2, 0, 1, true, true), 0, 1, false);
        lu = R::dnorm( R::qnorm(ld + log( 1 + exp ( lpu - ld) ) + log2, 0, 1, true, true), 0, 1, false );
        zs[i] = (lu - ll) / 2 / exp(ld);
        //Rprintf("%lg %lg %lg\n", ll, lu, ld);
      }
    }
    else { // Q_l < 0.5 < Q_u
      // (-dnorm(qnorm(2*ql))-dnorm(qnorm(2*(1-qu))))/2/(qu-ql) 
      ll = R::dnorm( R::qnorm(lpl + log2, 0, 1, true, true), 0, 1, false);
      lu = R::dnorm( R::qnorm(lpu + log2, 0, 1, true, true), 0, 1, false);
      zs[i] = (-lu - ll) / 2 / exp(ld);
    }
  }
  return(zs);
}
