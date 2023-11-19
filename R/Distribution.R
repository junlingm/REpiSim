#' An R6 class that defines a distribution to be used by a Calibrator object
#' 
#' @details Each Distribution object must implement a likelihood function to
#' be used by maximum likelihood calibrator to define the likelihood function.
#' A prior function must also be defined to provide a prior distirbution for
#' MCMC (such as used by the metrop method in the mcmc package)
Distribution <- R6Class(
  "Distribution",
  
  public = list(
    #' Provides the information to contruct a likelihood function.
    #' @return a named list with two components.
    #' 
    #' @details For the returned list, $par contains distribution 
    #' parameters that need to be specified in addition to the mean, 
    #' $logL contains a function that returns the log-likelihood, in the form of
    #' 
    #' $logL = function(x, mean, ...)
    #' 
    #' where ... contains necessary distribution parameters
    #' 
    likelihood = function() {
      stop("must be implemented by subclasses")
    },
    
    #' Provides a function to compute the log probability 
    #' @return a function with the formal function(x), which returns the log prior
    #' probability for a parameter value x
    prior = function(...) {
      stop("must be implemented by subclasses")
    }
  )
)

#' The Poisson distribution. 
#' 
#' @details This distribution has no extra parameter besides the mean
Poisson <- R6Class(
  "Poisson",
  inherit = Distribution,
  
  public = list(
    likelihood = function() {
      list(
        par = NULL,
        logL = function(x, mean, ...) {
          if (any(mean<0) || any(!is.finite(mean))) -Inf else
            sum(dpois(x, mean, log=TRUE))
        }
      )
    },
    
    prior = function(lambda) {
      function(x) dpois(x, lambda)
    }
  )
)

#' The negative binomial distribution. 
#' 
#' @details This distribution has one parameter named "size" in addition to the
#' mean, which is the size parameter in stats::dnbinom
NBinom <- R6Class(
  "NBinom",
  inherit = Distribution,
  
  public = list(
    likelihood = function() {
      list(
        par = "size",
        logL = function(x, mean, size) {
          if (any(mean<0) || any(!is.finite(mean)) || size<=0) -Inf else {
            sum(dnbinom(x, size=size, mu=mean, log=TRUE))
          }
        }
      )
    },
    
    #' Provides a function to compute the log probability 
    #' @param size the size parameter in stats::dnbinom
    #' @param prob the prob parameter in stats::dnbinom
    #' @return a function with the formal function(x), which returns the log prior
    #' probability for a parameter value x
    prior = function(size, prob) {
      function(x) dnbinom(x, shape, prob)
    }
  )
)

#' The normal distribution. 
#' 
#' @details This distribution has one parameter named "size" in addition to the
#' mean, which is the size parameter in stats::dnbinom
Normal <- R6Class(
  "Normal",
  inherit = Distribution,
  
  private = list(
    sd = NULL
  ),
  
  public = list(
    #' constructor
    #' @param sd the standard deviation
    #' @details 
    #' If sd is a number, it is used as a fixed standard deviation, the distribution
    #' has no parameter.
    #' 
    #' If sd is "fixed", the standard deviation is fixed, which is specified by
    #' the distribution parameter sd.
    #' 
    #' If sd is "proportional", the standard deviation is proportional to the mean,
    #' the coefficient is specified by a distribution parameter coef
    initialize = function(sd=c("fixed", "proportional")) {
      if (is.numeric(sd)) {
        private$sd = abs(sd)
      } else {
        private$sd = if (missing(sd)) "fixed" else
          match.arg(sd[[1]], c("fixed", "proportional"))
      }
    },
    
    likelihood = function() {
      if (is.numeric(private$sd)) {
        list(
          par = NULL,
          logL = function(x, mean) {
            if (any(!is.finite(mean))) -Inf else {
              sum(dnorm(x, mean=mean, sd=private$sd, log=TRUE))
            }
          }
        )
      } else {
        if (private$sd == "fixed") {
          list(
            par = "sd",
            logL = function(x, mean, sd) {
              if (any(!is.finite(mean)) || sd <=0 ) -Inf else {
                sum(dnorm(x, mean=mean, sd=sd, log=TRUE))
              }
            }
          )
        } else {
          list(
            par = "coef",
            logL = function(x, mean, coef) {
              if (any(!is.finite(mean)) || coef <= 0) -Inf else {
                sum(dnorm(x, mean=mean, sd=coef*abs(mean), log=TRUE))
              }
            }
          )
        }
      }
    },
    
    #' Provides a function to compute the log probability 
    #' @param mean the mean of the normal distribution
    #' @param sd the standard deviation of the normal distribution
    #' @return a function with the formal function(x), which returns the log prior
    #' probability for a parameter value x
    prior = function(mean, sd) {
      if (missing(sd)) {
        if (!is.numeric(private$sd)) stop("sd must be a number")
        sd = private$sd
      }
      function(x) dnorm(x, mean, sd)
    }
  )
)
