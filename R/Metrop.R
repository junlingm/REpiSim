#' A subclass of Baysian calibrator that uses the mcmc::metrop function
Metrop <- R6::R6Class(
  "Metrop",
  inherit = Baysian,
  
  private = list(
    .likelihood = NULL,
    .result = NULL,
    
    logL = function(pars, initial.values, parms, names) {
      n = names(private$.priors)
      names(pars) = names
      likelihood(private$.likelihood, private$.data, private$simulate, 
                 pars, initial.values, parms) +
        sum(sapply(1:length(n), function(i) {
          x = pars[[n[[i]]]]
          private$.priors[[i]]$log.density(x, pars)
        }))
    },
    
    fit = function(guess, initial.values, parms, ...) {
      x = metrop(private$logL, guess, ..., 
                 initial.values=initial.values, 
                 parms=parms, 
                 names=names(guess))
      private$.result = list(x=x, initial.values=initial.values, parms=parms)
      x
    },
    
    interpret = function(result) {
      result
    }
  ),
  
  public = list (
    #' @description initializer
    #' @param model the model to calibrate
    #' @param time either a numeric vector containing the times (including the 
    #' initial time) of the ODE solution that corresponds to the data, or a 
    #' character value giving the name of the column in data that corresponds 
    #' to time.
    #' @param data a data.frame object containign the data for the calibration
    #' @param likelihood a Distribution object specify the type of likelihood
    #' @param ... each argument is a formula defining the maps between 
    #' the data columns and the model variables. Please see the details section.
    #' @param cumulative whether the data is cumulative
    #' @details 
    #' A mapping is a named argument, where name is the
    #' data colummn name, and value corresponds to the model variables (or an 
    #' expression to calculate from the model variables.)
    initialize = function(model, time, data, likelihood, ..., cumulative=FALSE, mapping=character()) {
      library(mcmc)
      if (!"Likelihood" %in% class(likelihood)) 
        stop("likelihood must be a Likelihood object")
      private$.likelihood = likelihood
      super$initialize(model, time, data, ..., cumulative = cumulative, mapping = mapping)
    },
    
    continue = function(...) {
      if (is.null(private$.result))
        stop("was not run before")
      x = metrop(private$.result$x, ..., 
                 initial.values=private$.result$initial.values, 
                 parms=private$.result$parms,
                 names=names(private$.result$x$initial))
      private$.result$x = x
      private$interpret(x)
    }
  ),
  
  active = list(
    #' @field the names of all parameters
    parameters = function() {
      c(private$.model$parameters, private$.likelihood$par)
    }
  )
)
