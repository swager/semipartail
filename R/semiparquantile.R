#' Semiparametric estimate of quantiles of the main population
#'
#' @param main.sample the sample of interest
#' @param background.sample a (usually large) background sample used to stabilize tail inference
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param threshold threshold defining the beginning of the tail (selected automatically if not specified)
#'
#' @return plug-in estimates of the quantiles
semipar.quantile <- function(main.sample,
                         background.sample,
                         probs,
                         threshold=NULL) {
  # precompute semiparametric estimate of the cumulative distribution function
  semipar.est <- semipar.tail(main.sample, background.sample, threshold)
  cdf <- cumsum(semipar.est$weights)
  
  quantiles <- numeric(length(probs))
  for(i in 1:length(probs)) {
    p <- probs[i]
    if(p < 0 || p > 1) {
      stop("probs must be between 0 and 1")
    }
    
    # if p is in between two atoms, interpolate them
    if(p %in% cdf[1:(length(cdf)-1)]) {
      which.atom <- min(which(cdf == p))
      quantiles[i] <- mean(semipar.est$X[which.atom + 0:1])
    } 
    # otherwise, use the value for the first atom whose cdf exceeds p
    else {
      which.atom <- 1 + sum(cumsum(semipar.est$weights) < p)
      quantiles[i] <- semipar.est$X[which.atom]      
    }
  }
  return(quantiles)
}
