#
# Utility functions
#

#' Hill estimate of the tail index gamma as a function of threshold t
#'
#' @param x the raw data
#' @param t the threshold above which to compute the estimator
#'
#' @return a tail index estimate
hill <- function(x, t) {
  if (t <= 0 || sum(x >= t) <= 2) {
    return(NaN)
  }
  x.big = x[x >= t]
  n.big = length(x.big)
  gamma.hat = (sum(log(x.big)) - log(min(x.big))) / (n.big - 1)
  return (gamma.hat)
}

#' TODO

pareto.neg.loglik <- function(log.sigma, y, gma) {
    length(y)*log.sigma +
        (1+1/gma) * sapply(log.sigma, function(ls) sum(log(1+gma*y/exp(ls))))
}

gpd.scale <- function (x, t, gamma) {
  exp(optimize(function(ls) pareto.neg.loglik(ls, x[x > t], gamma), 
    interval=c(-20,20))$minimum)
}



#' Semiparametric estimate of the main.sample
#'
#' @param main.sample the sample of interest
#' @param background.sample a (usually large) background sample used to stabilize tail inference
#' @param threshold threshold defining the beginning of the tail (selected automatically if not specified)
#'
#' @return an estimated distribution, formatted as pairs
#'         (X = sample location, weights = amount of probability mass at X),
#'         where the weights sum to 1
semipar.tail <- function(main.sample,
                         background.sample,
                         threshold=NULL) {
              
    # If t is not passed in, use Guillou-Hall method to select threshold
    if(is.null(threshold)) {
       stop("Not implemented")
    }

    # If no big observations, return mean
    if(sum(main.sample > threshold)==0) {
      return(data.frame(X = main.sample,
             weights = 1/length(main.sample)))
    }

    # Step 1: Estimate gamma and sigma for the background distribution
    # to get scale of sufficient statistic T
    
    gamma.hat <- hill(background.sample, threshold)
    sigma.hat <- gpd.scale(background.sample, threshold, gamma.hat)
    kappa <- sigma.hat/gamma.hat

    # Step 2: Estimate eta (logistic regression slope)
    
    # Semiparametric inference only uses the observations above the threshold
    big.x <- main.sample[main.sample > threshold] - threshold
    big.y <- background.sample[background.sample > threshold] - threshold
    
    # Compute the logistic regression features
    Tx <- big.x / (big.x + kappa)
    Ty <- big.y / (big.y + kappa)
    
    # Estimate eta via logistic regression
    Tz <- c(Tx, Ty)    
    response <- rep(1:0, c(length(Tx), length(Ty)))    
    eta.hat <- coef(glm(response ~ Tz, family = binomial()))[2]

    # Step 3: Estimate the full cdf
    
    # below threshold, use empirical cdf
    low.data <- data.frame(X = main.sample[main.sample <= threshold],
                           weights = 1 / length(main.sample))
                 
    # above threshold, tilt background cdf
    background.rel.weights <- exp(eta.hat * Ty) / sum(exp(eta.hat * Ty))
    background.weights <- background.rel.weights * sum(main.sample > threshold) / length(main.sample)
    high.data <- data.frame(X = background.sample[background.sample > threshold],
                            weights = background.weights)
    
    # combine the two distribution estimates into a single final answer
    semipar.sample <- rbind(low.data, high.data)
    return(semipar.sample)
}
