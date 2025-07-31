#' An R6 class that defines a likelihood function to be used by a Calibrator object
#' 
#' @details Each likelihood object must provide a read-only active field
#' name par specifying likelihood parameters, and function logL for log
#' likelihood. These must be implemented by subclasses
#' @name Likelihood
#' @docType class
#' @export
#' 

Likelihood <- R6Class(
  "Likelihood",
  public = list(
    #' @description The log likelihood function
    #' @param x the observation to calculate the likelihood
    #' @param mean the mean of the likelihood distribution
    #' @param ... other likelihood parameter values
    #' @return a numeric value giving the likelihood.
    #' 
    #' @details x and mean can be vectors
    logL = function(x, mean, ...) {
      stop("must be implemented by subclasses")
    }
  ),
  
  active = list(
    #' @field par return a vector of parameter names for this likelihood distribution
    par = function() {
      NULL
    }
  )
)

#' The Poisson distribution
#' 
#' @details This distribution has no extra parameter besides the mean
#' @name Poisson
#' @docType class
#' @export
#' 
Poisson <- R6Class(
  "Poisson",
  inherit = Likelihood,
  
  public = list(
    #' @description The log likelihood function for the Poisson distribution.
    #' @param x the value to calculate the likelihood
    #' @param mean the mean of the Poisson distribution
    #' @param ... other parameters (not used)
    logL = function(x, mean, ...) {
      if (any(mean<0) || any(!is.finite(mean))) -Inf else
        sum(dpois(x, mean, log=TRUE))
    }
  )
)

#' The negative binomial distribution
#' 
#' @details This distribution has one parameter named "size" in addition to the
#' mean, which is the size parameter in stats::dnbinom
#' @name NBinom
#' @docType class
#' @export
#' 

NBinom <- R6Class(
  "NBinom",
  inherit = Likelihood,
  
  public = list(
    #' @description The log likelihood function for the negative binomial distribution
    #' @param x the value to calculate the likelihood
    #' @param mean the mean of the negative binomial distribution
    #' @param size the size parameter of the negative binomial distribution
    logL = function(x, mean, prob) {
      if (any(mean<0) || any(!is.finite(mean)) || prob<=0) -Inf else {
        size = mean * prob / (1-prob)
        sum(dnbinom(x, size=size, prob=prob, log=TRUE))
      }
    }
  ),
  
  active = list(
    #' @field par return a vector of parameter names for this likelihood
    par = function() {
      "prob"
    }
  )
)

#' The normal distribution
#' 
#' @name Normal
#' @docType class
#' @export
#' 

Normal <- R6Class(
  "Normal",
  inherit = Likelihood,
  
  public = list(
    #' @description The log likelihood function for the normal distribution
    #' @param x the value to calculate the likelihood
    #' @param mean the mean of the normal distribution
    #' @param sd the standard deviation of the normal distribution
    logL = function(x, mean, sd) {
      if (any(!is.finite(mean)) || sd <=0 ) -Inf else {
        sum(dnorm(x, mean=mean, sd=sd, log=TRUE))
      }
    }
  ),
  
  active = list(
    #' @field par return a vector of parameter names for this likelihood
    par = function() {
      "sd"
    }
  )
)


#' A normal distribution which standard deviation is proportion to mean
#' 
#' @details This distribution has one parameter named "ceof" in addition to the
#' mean, where sd = coef*mean
#' @name NormalProportional
#' @docType class
#' @export
#' 

NormalProportional <- R6Class(
  "NormalProportional",
  inherit = Likelihood,
  
  public = list(
    #' @title The log likelihood function for the normal distribution
    #' @param x the value to calculate the likelihood
    #' @param mean the mean of the normal distribution
    #' @param coef the coefficient of variation, where sd = coef*mean
    logL = function(x, mean, coef) {
      if (any(!is.finite(mean)) || coef <= 0) -Inf else {
        sum(dnorm(x, mean=mean, sd=coef*abs(mean), log=TRUE))
      }
    }
  ),
  
  active = list(
    #' @field par return a vector of parameter names for this likelihood
    par = function() {
      "coef"
    }
  )
)
