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


#' expectation-based beta-binomial overdispersion test
#
#' @param xs A vector of non-negative counts representing observed data
#' @param sizes A vector of positive values representing total counts
#' @param ws A vector of weights for each observed data
#' @param alphas A non-negative vector of parameters of the beta distribution
#' @param betas A non-negative vector of parameters of the beta distribution
#' @param probs (alternative) A vector of binomial probabilities for each observed data
#' @param rhos (alternative) A vector of beta-binomial overdispersion parameters rho = 1/(alpha+beta+1)
#' @param pos.only Ignore zero observations
#' @param adj.var Adjust variance to improved power
#' @param approx.under Use approximate variance calculation when Pr(X)<value
#' @param cap.z The threshold that an individual z-score can contribute to the test statistics
#' @param var.add When adj.var=TRUE, amount of variance to add to the denominator to prevent anti-conservative behavior due to variance adjustment
#' @return z-score from the meta-analysis
#' @export
betabinom_overdisp_test_r_n <- function(xs, sizes,ws=NULL,alphas=NULL, betas=NULL,
                                      probs = NULL, rhos = NULL, 
                                      pos.only=TRUE, adj.var=TRUE, 
                                      approx.under = 1e-4, cap.z=10, var.add=1) {
  if (is.null(ws) | length(ws)==1) {
    ws = rep(1, length(xs))
  }
  if ( pos.only ) {
    iv <- (xs > 0)
    xs <- xs[iv]
    sizes <- sizes[iv]
    alphas <- alphas[iv]
    betas <- betas[iv]
    ws <- ws[iv]
  }
  if ( adj.var == FALSE ) {
    var.add <- 0
  }
  n = length(xs)
  fadj = rep(0, n)
  if ( pos.only ) {
    fadj = cpp_pbbinom(rep(0,n),sizes,alphas,betas,0,1)
  }
  if ( is.null(alphas) | is.null(betas) ) {
    if ( is.null(probs) | is.null(rhos) ) {
      stop("At least one set of beta binomial parameters, (alpha,beta) or (probs, rho), has to be provided.")
    }
    alphas = probs * (1-rhos)/rhos
    betas  = (1-probs) * (1-rhos)/rhos
  }
  
  qls = cpp_pbbinom(xs,   sizes, alphas, betas, 0, 1) - fadj
  qus = cpp_pbbinom(xs-1, sizes, alphas, betas, 0, 1) - fadj
  zs = multi_bidir_zs_from_qt_n(qls, qus, 1, 1)
  
  mids = round(sizes * alphas / (betas + alphas))
  minx = rep(as.integer(pos.only),n)
  maxx = sizes
  if (approx.under > 0) {
    for (i in 1:n) {
      for (j in (mids[i]+1):n) {
        if (cpp_dbbinom(j,sizes[i],alphas[i],betas[i],0) < approx.under) {
          maxx[i] = j
          break
        }
      }
      for (j in (mids[i]:0)) {
        if (cpp_dbbinom(j,sizes[i],alphas[i],betas[i],0) < approx.under) {
          minx[i] = j
          break
        }
      }
    }    
    # minx[minx < pos.only] = pos.only
  } 
  vs = sapply(1:n, function(i) {
    qt = cpp_pbbinom(maxx[i]:minx[i], sizes[i], alphas[i], betas[i], 0, 1)
    qt = qt - fadj[i]
    if (qt[length(qt)] < log(0.5) ) {
      qt = c(qt, 0)
    }    
    return(bidir_etn_var_from_qt( qt, 1, 1 ))
  })
  
  zs[zs > cap.z] = cap.z
  zs[zs < 0-cap.z] = 0-cap.z
  return( sum(ws * zs) / sqrt(var.add * mean(ws^2) + sum(vs * ws^2)) )
}



#' Expectation-based beta-binomial overdispersion test
#
#' @param xs A vector of non-negative counts representing observed data
#' @param sizes A vector of positive values representing total counts
#' @param ws A vector of weights for each observed data
#' @param alphas A non-negative vector of parameters of the beta distribution
#' @param betas A non-negative vector of parameters of the beta distribution
#' @param probs (alternative) A vector of binomial probabilities for each observed data
#' @param rhos (alternative) A vector of beta-binomial overdispersion parameters rho = 1/(alpha+beta+1)
#' @param pos.only Ignore zero observations
#' @param adj.var Adjust variance to improved power
#' @param approx.under Use approximate variance calculation when Pr(X)<value
#' @param cap.z The threshold that an individual z-score can contribute to the test statistics
#' @param var.add When adj.var=TRUE, amount of variance to add to the denominator to prevent anti-conservative behavior due to variance adjustment
#' @return z-score from the meta-analysis
#' @export
betabinom_overdisp_test_n <- function(xs, sizes,ws=NULL,alphas=NULL, betas=NULL,
                                      probs = NULL, rhos = NULL,
                                      pos.only=TRUE, adj.var=TRUE,
                                      approx.under = 1e-4, cap.z=10, var.add=1) {
  if (is.null(ws) | length(ws)==1) {
    ws = rep(1, length(xs))
  }
  if ( pos.only ) {
    iv <- (xs > 0)
    xs <- xs[iv]
    sizes <- sizes[iv]
    alphas <- alphas[iv]
    betas <- betas[iv]
    ws <- ws[iv]
  }
  if ( adj.var == FALSE ) {
    var.add <- 0
  }
  n = length(xs)
  df <- suppressWarnings(betabinom_multi_bidir_n(xs, sizes, alphas, betas,
                                                 pos_only=pos.only, var_adj=adj.var,
                                                 approx_under=approx.under))
  df$zs[df$zs > cap.z] <- cap.z
  df$zs[df$zs < -cap.z] <- -cap.z
  return( sum(ws * df$zs) / sqrt(var.add * mean(ws^2) + sum(df$variance * ws^2)) )
}

