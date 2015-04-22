#' Semiparametric estimate of the expectation of the main population
#'
#' @param main.sample the sample of interest
#' @param background.sample a (usually large) background sample used to stabilize tail inference
#'
#' @return numeric estimate of the expectation
semipar.mean <- function(main.sample,
                         background.sample,
                         threshold=NULL) {
  with(semipar.tail(main.sample, background.sample, threshold),
       sum(X * weights) / sum(weights))
}
