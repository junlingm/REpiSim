#' A calibrator that uses a Baysian method.
#' @name Baysian
#' @docType class
#' @export
#' 

Baysian <- R6::R6Class(
  "Baysian",
  inherit = Calibrator,
  
  private = list(
    .priors = NULL,
    
    simulator = function(model) {
      ODE$new(model)
    },
    
    fit = function(guess, initial.values, parms, ...) {
      stop("must be implemented by subclasses")
    },
    
    prior = function(formula) {
      Prior$new(formula)
    },
    
    .calibrate = function(pars, initial.values, parms, priors, guess, ...) {
      if (!is.list(priors))
        stop("invalid priors")
      np = names(priors)
      if (is.null(np) || any(np=="")) 
        stop("priors must be named")
      extra = setdiff(np, pars)
      if (length(extra) > 0)
        stop("extra prior", if(length(extra)==1) "" else "s", 
             ": ", paste(extra, collapse=", "))
      miss = setdiff(pars, np)
      if (length(miss) > 0)
        stop("missing prior", if(length(miss)==1) "" else "es", ": ", 
             paste(miss, collapse=", "))
      private$.priors = lapply(priors, private$prior)
      private$fit(guess, initial.values=initial.values, parms=parms, ...)
    }
  ),
  
  public = list (
    #' @title Calibrate the model to data
    #' 
    #' @param initial.values the initial values for the model. The parameters 
    #' that need to be estimate should be NA, those that do not need to be
    #' estimated must contain a finite value.
    #' @param parms the parameter values of the model. The parameters 
    #' that need to be estimate should be NA, those that do not need to be
    #' estimated must contain a finite value.
    #' @param priors a named list of Distribution objects specifying the priors.
    #' Each name is a parameter.
    #' @param guess the initial guess of the parameters to be fitted
    #' @param ... extra arguments to be passed to calibrators
    calibrate = function(initial.values, parms, priors, guess, ...) {
      super$calibrate(initial.values, parms, priors, guess, ...)
    }
  ),
  
  active = list(
    #' @field samples the samples from the fit, a read-only field
    samples = function() {
      NULL
    }
  )
)
