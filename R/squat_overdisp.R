#' binom_overdisp_test_n
#' expectation-based binomial overdispersion test
#
#' @param xs A vector of non-negative counts representing observed data
#' @param sizes A vector of positive values representing total counts
#' @param ps A vector of binomial probabilities for each observed data
#' @param ws A vector of weights for each observed data
#' @param pos.only Ignore zero observations
#' @param adj.var Adjust variance to improved power
#' @param approx.under Use approximate variance calculation when Pr(X)<value
#' @param cap.z The threshold that an individual z-score can contribute to the test statistics
#' @param var.add When adj.var=TRUE, amount of variance to add to the denominator to prevent anti-conservative behavior due to variance adjustment
#' @return z-score from the meta-analysis
#' @export
binom_overdisp_test_n <- function(xs, sizes, ps, ws, pos.only=TRUE, adj.var=TRUE, approx.under = 1e-4, cap.z=10, var.add=1) {
  if ( pos.only ) {
    iv <- (xs > 0)
    xs <- xs[iv]
    sizes <- sizes[iv]
    ps <- ps[iv]
    ws <- ws[iv]
  }
  if ( adj.var == FALSE ) {
    var.add <- 0
  }
  df <- suppressWarnings(binom_multi_bidir_n(xs, sizes, ps, pos.only, adj.var, approx.under))
  df$zs[df$zs > cap.z] <- cap.z
  df$zs[df$zs < -cap.z] <- -cap.z
  return( sum(ws * df$zs) / sqrt(var.add * mean(ws^2) + sum(df$variance * ws^2)) )
}

#' betabinom_overdisp_test_n
#' expectation-based beta-binomial overdispersion test
#
#' @param xs A vector of non-negative counts representing observed data
#' @param sizes A vector of positive values representing total counts
#' @param probs A vector of binomial probabilities for each observed data
#' @param obs A vector of beta-binomial overdispersion parameters = (1-rho)/rho = alpha+beta
#' @param ws A vector of weights for each observed data
#' @param pos.only Ignore zero observations
#' @param adj.var Adjust variance to improved power
#' @param approx.under Use approximate variance calculation when Pr(X)<value
#' @param cap.z The threshold that an individual z-score can contribute to the test statistics
#' @param var.add When adj.var=TRUE, amount of variance to add to the denominator to prevent anti-conservative behavior due to variance adjustment
#' @return z-score from the meta-analysis
#' @export
betabinom_overdisp_test_n <- function(xs, sizes, probs, ovs, ws, 
                                      pos.only=TRUE, adj.var=TRUE, 
                                      approx.under = 1e-4, cap.z=10, var.add=1) {
  if ( pos.only ) {
    iv <- (xs > 0)
    xs <- xs[iv]
    sizes <- sizes[iv]
    probs <- probs[iv]
    ovs <- ovs[iv]
    ws <- ws[iv]
  }
  if ( adj.var == FALSE ) {
    var.add <- 0
  }
  n = length(xs)
  fadj = rep(1, n)
  if ( pos.only ) {
    fadj = 1-dbetabinom(rep(0,n),sizes,probs,ovs)
  }
  qls = (1-pbetabinom(xs,   sizes, probs, ovs))/fadj
  qus = (1-pbetabinom(xs-1, sizes, probs, ovs))/fadj
  zs = multi_bidir_zs_from_qt_n(qls, qus, 0, 1)
  minx = rep(as.integer(pos.only),n)
  maxx = sizes
  if (approx.under > 0) {
    minx = qbetabinom(approx.under, sizes, probs, ovs);
    maxx = qbetabinom(1-approx.under, sizes, probs, ovs);
    minx[minx < pos.only] = pos.only
  } 
  vs = sapply(1:n, function(i) {
    if (minx[i] == maxx[i]) {
      qt = 1-pbetabinom(maxx[i],   sizes[i], probs[i], ovs[i])
    } else {
      qt = 1-pbetabinom(maxx[i]:minx[i], sizes[i], probs[i], ovs[i])
    }
    qt = qt/fadj[i]
    if (qt[length(qt)] < 0.5) {
      qt = c(qt, 1)
    }    
    return(bidir_etn_var_from_qt( qt, 0, 1 ))
  })
  zs[zs > cap.z] = cap.z
  zs[zs < 0-cap.z] = 0-cap.z
  return( sum(ws * zs) / sqrt(var.add * mean(ws^2) + sum(vs * ws^2)) )
}