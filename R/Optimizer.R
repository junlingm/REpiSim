#' A subclass of Calibrator that uses an optimizer and an initial guess.
Optimizer <- R6::R6Class(
  "Optimizer",
  inherit = Calibrator,
  
  private = list(
    simulator = function(model) {
      ODE$new(model)
    },
    
    objective = function(pars, ...) {
      stop("The objective method must be implemented by a subclass")
    },
    
    optimizer = function(pars, initial.values, parms, ...) {
      stop("The optimizer method must be implemented by a subclass")
    },
    
    .calibrate = function(pars, initial.values, parms, guess, ...) {
      if (!is.numeric(guess) || any(is.na(guess)))
        stop("invalid initial values")
      ng = names(guess)
      if (is.null(ng) || any(ng=="")) 
        stop("initial guesses must be named")
      extra = setdiff(pars, ng)
      if (length(extra) > 0)
        stop("variable", if(length(extra)==1) "" else "s", 
             " not defined in model: ", paste(extra, collapse=", "))
      miss = setdiff(ng, pars)
      if (length(miss) > 0)
        stop("missing initial guess", if(length(miss)==1) "" else "es", ": ", 
             paste(miss, collapse=", "))
      private$optimizer(guess[pars], initial.values=initial.values, parms=parms, ...)
    }
  ),
  
  public = list (
    #' Calibrate the model to data
    #' 
    #' @param initial.values the initial values for the model. The parameters 
    #' that need to be estimate should be NA, those that do not need to be
    #' estimated must contain a finite value.
    #' @param parms the parameter values of the model. The parameters 
    #' that need to be estimate should be NA, those that do not need to be
    #' estimated must contain a finite value.
    #' @param guess the initial guess of the parameters to be fitted
    #' @param ... extra arguments to be passed to calibrators
    calibrate = function(initial.values, parms, guess, ...) {
      super$calibrate(initial.values, parms, guess, ...)
    }
  )
)
