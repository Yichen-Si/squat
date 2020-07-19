#' squat_binom_overdisp_rand : randomized binomial overdispersion test
#
#' @param xs A vector of non-negative counts representing observed data
#' @param sizes A vector of positive values representing total counts
#' @param ps A vector of binomial probabilities for each observed data
#' @param ws A vector of weights for each observed data
#' @param nrep Number of repetitions
#' @param pos.only Ignore zero observations
#' @return A (nrep x 2) matrix, containing the following values in each column
#' * Z-score from randomzied exact meta-analysis
#' * -log10 p-value corresponding to the z-score
#' @export
squat_binom_overdisp_rand <- function(xs, sizes, ps, ws, nrep=1, pos.only=TRUE) {
  if ( pos.only ) {
    iv <- (xs > 0)
    xs <- xs[iv]
    sizes <- sizes[iv]
    ps <- ps[iv]
    ws <- ws[iv]
    k <- length(xs)
    f1ps <- 1 - (1-ps)^sizes
  }
  else {
    k <- length(xs)
    f1ps <- rep(1,length(xs))
  }
  dx <- dbinom(xs, sizes, ps)  # dx represent binomial density in log-scale
  px <- pbinom(xs, sizes, ps, lower.tail=FALSE) # px represents binomial upper-tail cdf
  rnd.mtx <- ( matrix(px, k, nrep) + matrix(dx, k, nrep) * matrix(runif(k*nrep), k, nrep) ) / matrix(f1ps, k, nrep) * 2
  idx.flip <- (rnd.mtx > 1)
  rnd.mtx[idx.flip] <- 2 - rnd.mtx[idx.flip]
  rnd.mtx[rnd.mtx < 1e-30] <- 1e-30
  zs <- colSums(matrix(-ws, k, nrep) * qnorm(rnd.mtx) / sqrt(sum(ws^2)))
  logps <- 0-pnorm(zs, lower.tail=FALSE,log.p=TRUE)/log(10)
  return(cbind(zs,logps))
}

#' squat_binom_overdisp_test : expectation-based binomial overdispersion test
#
#' @param xs A vector of non-negative counts representing observed data
#' @param sizes A vector of positive values representing total counts
#' @param ps A vector of binomial probabilities for each observed data
#' @param ws A vector of weights for each observed data
#' @param pos.only Ignore zero observations
#' @param adj.var Adjust variance to improved power
#' @param approx.under Use approximate variance calculation when Pr(X)<value
#' @param cap.z The threshold that an individual z-score can contribute to the test statistics
#' @return z-score from the meta-analysis
#' @export
squat_binom_overdisp_test <- function(xs, sizes, ps, ws, pos.only=TRUE, adj.var=TRUE, approx.under = 1e-4, cap.z=10) {
  if ( pos.only ) {
    iv <- (xs > 0)
    xs <- xs[iv]
    sizes <- sizes[iv]
    ps <- ps[iv]
    ws <- ws[iv]
  }
  zs <- suppressWarnings(squat_multi_binom_bidir(xs, sizes, ps, pos.only, adj.var, approx.under))
  zs[zs > cap.z] <- cap.z
  zs[zs < -cap.z] <- -cap.z
  return( sum(ws * zs) / sqrt(sum(ws^2)) )
}

#' squat_binom_directional_test : expectation-based binomial directional test
#
#' @param xs A vector of non-negative counts representing observed data
#' @param ys A vector of 0/1 (or -1/1) representing the direction of test
#' @param sizes A vector of positive values representing total counts
#' @param ps A vector of binomial probabilities for each observed data
#' @param ws A vector of weights for each observed data
#' @param adj.var Adjust variance to improved power
#' @param approx.under Use approximate variance calculation when Pr(X)<value
#' @param cap.z The threshold that an individual z-score can contribute to the test statistics
#' @return z-score from the meta-analysis
#' @export
squat_binom_directional_test <- function(xs, ys, sizes, ps, ws, adj.var=TRUE, approx.under = 1e-4, cap.z = 10) {
  zs <- suppressWarnings(squat_multi_binom_unidir(xs, sizes, ps, adj.var, approx.under))
  zs[zs > cap.z] <- cap.z
  zs[zs < -cap.z] <- -cap.z
  ys[ys == 0] <- -1
  return ( sum(ws * ys * zs) / sqrt(sum(ws^2)) )
}
