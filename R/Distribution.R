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
