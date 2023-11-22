#' An R6 class that defines a likelihood function to be used by a Calibrator object
#' 
#' @details Each likelihood object must provide a read-only active field
#' name par specifying likelihood parameters, and function logL for log
#' likelihood. These must be implemented by subclasses
Likelihood <- R6Class(
  "Likelihood",
  
  public = list(
    #' The log likelihood function.
    #' 
    #' @param x the observation to calculate the likelihood
    #' @param mean the mean of the likelihood distribution
    #' @param ... other likelihood parameter values
    #' @return a numeric value giving the likelihood.
    #' 
    #' @detail x and mean can be vectors
    logL = function(x, mean, ...) {
      stop("must be implemented by subclasses")
    }
  ),
  
  active = list(
    #' @field par return a vector of parameter names for this likelihood 
    #' distribution
    par = function() {
      NULL
    }
  )
)

#' The Poisson distribution. 
#' 
#' @details This distribution has no extra parameter besides the mean
Poisson <- R6Class(
  "Poisson",
  inherit = Likelihood,
  
  public = list(
    logL = function(x, mean, ...) {
      if (any(mean<0) || any(!is.finite(mean))) -Inf else
        sum(dpois(x, mean, log=TRUE))
    }
  )
)

#' The negative binomial distribution. 
#' 
#' @details This distribution has one parameter named "size" in addition to the
#' mean, which is the size parameter in stats::dnbinom
NBinom <- R6Class(
  "NBinom",
  inherit = Likelihood,
  
  public = list(
    logL = function(x, mean, size) {
      if (any(mean<0) || any(!is.finite(mean)) || size<=0) -Inf else {
        sum(dnbinom(x, size=size, mu=mean, log=TRUE))
      }
    }
  ),
  
  active = list(
    par = function() {
      "size"
    }
  )
)

#' The normal distribution. 
#' 
#' @details This distribution has one parameter named "size" in addition to the
#' mean, which is the size parameter in stats::dnbinom
Normal <- R6Class(
  "Normal",
  inherit = Likelihood,
  
  public = list(
    logL = function(x, mean, sd) {
      if (any(!is.finite(mean)) || sd <=0 ) -Inf else {
        sum(dnorm(x, mean=mean, sd=sd, log=TRUE))
      }
    }
  ),
  
  active = list(
    par = function() {
      "sd"
    }
  )
)

#' A normal distribution which standard deviation is proportion to mean
#' 
#' @details This distribution has one parameter named "ceof" in addition to the
#' mean, where sd = coef*mean
NormalProportional <- R6Class(
  "NormalProportional",
  inherit = Likelihood,
  
  public = list(
    logL = function(x, mean, coef) {
      if (any(!is.finite(mean)) || coef <= 0) -Inf else {
        sum(dnorm(x, mean=mean, sd=coef*abs(mean), log=TRUE))
      }
    }
  ),
  
  active = list(
    par = function() {
      "coef"
    }
  )
)
