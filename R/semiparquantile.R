#' Semiparametric estimate of quantiles of the main population
#'
#' @param main.sample the sample of interest
#' @param background.sample a (usually large) background sample
#'        used to stabilize tail inference
#' @param quantiles numeric vector of probabilities with values in [0,1].
#' @param threshold threshold defining the beginning of the tail
#'        (selected automatically if not specified)
#'
#' @return Plug-in estimates of the specified quantiles for the
#'         main sample of interest.
#'
#' @export semipar.quantile
semipar.quantile <- function(main.sample,
                         background.sample,
                         quantiles,
                         threshold=NULL) {
  # precompute semiparametric estimate of the cumulative distribution function
  semipar.est <- semipar.tail(main.sample, background.sample, threshold)
  semipar.ord <- semipar.est[order(semipar.est$X),]
  semipar.cdf <- cumsum(semipar.ord$weights)
  
  quant.hat <- numeric(length(quantiles))
  for(i in 1:length(quantiles)) {
    p <- quantiles[i]
    if(p <= 0 || p >= 1) {
      stop("quantiles must be between 0 and 1")
    }
    if(p < min(semipar.cdf) || p >= max(semipar.cdf)) {
      stop("quantiles outside of range of estimated distribution")
    }
    
    which.atom <- 1 + sum(semipar.cdf < p)
    quant.hat[i] <- semipar.ord$X[which.atom]   
  }
  return(quant.hat)
}
