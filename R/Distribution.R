#' An R6 class that defines a distribution to be used as either a likelihood or a prior.
#' 
#' @details Each Distribution object must provide two methods, a "log.density" method
#' which is the log-density function. In addition, a "log.likelihood" method which is
#' the log-likelihood function and the distirbution parameters
#' @name Distribution
#' @docType class
#' @export

Distribution <- R6Class(
  "Distribution",
  private = list(
    ## the log likelihood
    .log.likelihood = NULL,
    ## the log density
    .log.density = NULL
  ),
  public = list(
    #' @param ... the parameters needed to specify the distribution. 
    #' @param .density the density function
    #' @details to use this object as a likelihood, additional parameters that are needed must
    #' be provided as NA, while the other parameters that needs to be inferred from the mean
    #' and the estimated parameters must be provided as a formula
    initialize = function(..., .density) {
      if (is.null(.density))
        stop(".density must be provided to initialize a Distribution object")
      args = list(...)
      ns = names(args)
      if (is.null(ns) || any(ns=="")) 
        stop("Distribution parameters must be named")
      # is this a distribution? if so, then the arguments in ... must all be numeric
      if (all(is.numeric(unlist(args)))) {
        call.args = c(x=as.name("x"), args, log=TRUE)
        private$.log.density = function(x) {
          do.call(.density, call.args)
        }
      } else { # otherwise, this is a likelihood
        call.args = list(x=NA, mean=NA)
        par = list()
        extra = character() # parameters in formula
        expr = list() #$ the list of formula defining parameters
        for (n in ns) {
          v = args[[n]]
          if (is.call(v)) {
            e = Expression$new(v)
            expr[[n]] = e
            extra = c(extra, e$parms)
            par[[n]] = as.name(n)
          } else if (is.name(v)) {
            extra = c(extra, v)
            par[[n]] = as.name(n)
          } else if (is.na(v)) {
            if (n != "mean") call.args[[n]] = NA
            par[[n]] = as.name(n)
          } else par[[n]] = v
        }
        expr = expr[order(expr)]
        # statements calculating parameters
        stmt = sapply(names(expr), function(n) {
          e = expr[[n]]
          call("<-", as.name(n), e$expr)
        })
        # add in parameters in the formula
        extra = setdiff(extra, ns)
        if (length(extra) > 0) call.args[extra] = NA
        density = as.call(c(
          list(as.name(".density"), as.name("x")), 
          par, log=TRUE
        ))
        sum = call("sum", density)
        body <- as.call(c(as.name("{"), stmt, sum))
        call.args = c(call.args, body)
        private$.log.likelihood = as.function(call.args)
      }
    }
  ),
  active = list(
    #' @field log.likelihood The log-likelihood function, which takes as least two parameters:
    #' 
    #' x the observation to calculate the likelihood
    #' 
    #' mean the mean of the distribution. This typically corresponds to the model solution
    #' 
    #' ... the parameters of the distribution,  in addition to the mean, to calculate the 
    #' likelihood
    #' 
    #' It return the log likelihood
    log.likelihood = function() {
      private$.log.likelihood
    },

    #' @field log.density the log-density function that takes a single paramter x,  which is
    #' the value where the density is evaluated at.
    #' it returns the log-density at x
    log.density = function() {
      private$.log.density
    }
  )
)

#' The Poisson distribution
#' 
#' @name Poisson
#' @docType class
#' @export

Poisson <- R6Class(
  "Poisson",
  inherit = Distribution,
  
  public = list(
    #' @param lambda the mean of the Poisson distribution. When getting the log-likelihood,
    #' this parameter is unused.
    initialize = function(lambda) {
      if (missing(lambda) || is.na(lambda)) lambda = quote(mean)
      super$initialize(lambda = lambda, .density = dpois)
    }
  )
)

#' The negative binomial distribution
#' 
#' @name NBinom
#' @docType class
#' @export

NBinom <- R6Class(
  "NBinom",
  inherit = Distribution,
  
  public = list(
    #' @param size the size parameter of the negative binomial distribution
    #' @param prob the probability parameter of the negative binomial distribution
    initialize = function(size, prob) {
      if (missing(size)) {
        size = quote(mean*prob/(1-prob))
      } else if (missing(prob)) {
        prob = quote(size/(mean+size))
      }
      super$initialize(size = size, prob = prob, .density = dnbinom)
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
  inherit = Distribution,
  
  public = list(
    #' @param mean the mean of the normal distribution
    #' @param sd the standard deviation of the normal distribution
    initialize = function(mean, sd) {
      if (missing(mean)) {
        mean = quote(mean)
      } else if (missing(sd)) {
        stop("the standard deviation must be provided for the normal distribution")
      }
      super$initialize(mean = mean, sd = sd, .density = dnorm)
    }
  )
)

#' The uniform distribution
#' @name Uniform
#' @docType class
#' @export

Uniform <- R6Class(
  "Uniform",
  inherit = Distribution,
  
  public = list(
    #' @param min the minimum value of the uniform distribution
    #' @param max the maximum value of the uniform distribution
    initialize = function(min, max) {
      if (missing(min)) {
        min = quote(2*mean - max)
      } else if (missing(max)) {
        max = quote(2*mean + min)
      }
      super$initialize(min = min, max = max, .density = dunif)
    }
  )
)

#' The exponential distribution
#' @name Exponential
#' @docType class
#' @export

Exponential <- R6Class(
  "Exponential",
  inherit = Distribution,
  
  public = list(
    #' @param rate the rate parameter of the exponential distribution
    initialize = function(rate) {
      if (missing(rate)) {
        rate = quote(1/mean)
      }
      super$initialize(rate = rate, .density = dexp)
    }
  )
)

#' The gamma distribution
#' @name Gamma
#' @docType class
#' @export

Gamma <- R6Class(
  "Gamma",
  inherit = Distribution,
  
  public = list(
    #' @param shape the shape parameter of the gamma distribution
    #' @param scale the scale parameter of the gamma distribution, i.e., 1/rate
    initialize = function(shape, scale) {
      if (missing(shape)) {
        shape = quote(mean/scale)
      } else if (missing(scale)) {
        scale = quote(mean/shape)
      }
      super$initialize(shape = shape, scale = scale, .density = dgamma)
    }
  )
)

#' The beta distribution
#' @name Beta
#' @docType class
#' @export

Beta <- R6Class(
  "Beta",
  inherit = Distribution,
  
  public = list(
    #' @param shape1 the first shape parameter of the beta distribution
    #' @param shape2 the second shape parameter of the beta distribution
    initialize = function(shape1, shape2) {
      if (missing(shape1)) {
        # mean = shape1/(shape1 + shape2). calculate shape1 in terms of mean and shape2
        shape1 = quote(mean/(1-mean)*shape2)
      } else if (missing(shape2)) {
        shape2 = quote((1-mean)/mean*a)
      }
      super$initialize(shape1 = shape1, shape2 = shape2, .density = dbeta)
    }
  )
)
