#' A subclass of Calibrator that uses an optimizer and an initial guess.
Optimizer <- R6::R6Class(
  "Optimizer",
  inherit = Calibrator,
  
  private = list(
    simulator = function(model) {
      ODE$new(model)
    },
    
    #' The objective function for optimization
    #' This method must be implemented by subclasses.
    #' @param par the values of the fitted parameters
    #' @param formula gives the parameter values that should be calculated
    #' using this formula
    #' @param fixed the parameter values that are given
    #' @return the fitting results, which will be passed to 
    #' the intepret method.
    objective = function(par, formula, fixed, ...) {
      stop("The objective method must be implemented by a subclass")
    },
    
    #' The optimization is done in this function.
    #' This method must be implemented by subclasses.
    #' @param guess the initial guesses of the fitted parameters
    #' @param formula gives the parameter values that should be calculated
    #' using this formula
    #' @param fixed the parameter values that are given
    #' @return the fitting results, which will be passed to 
    #' the intepret method.
    optimizer = function(guess, formula, fixed, ...) {
      stop("The optimizer method must be implemented by a subclass")
    },
    
    #' The actual calibration is done in this function.
    #' This method Must be implemented by subclasses.
    #' @param fit the parameter values that should be fitted
    #' @param formula gives the parameter values that should be calculated
    #' using this formula
    #' @param fixed the parameter values that are given
    #' @param guess the initial guess of the parameters to be fitted
    #' @return the fitting results, which will be passed to 
    #' the intepret method.
    .calibrate = function(fit, formula, fixed, guess, ...) {
      if (!is.numeric(guess) || any(is.na(guess)))
        stop("invalid initial values")
      ng = names(guess)
      if (is.null(ng) || any(ng=="")) 
        stop("initial guesses must be named")
      extra = setdiff(fit, ng)
      if (length(extra) > 0)
        stop("variable", if(length(extra)==1) "" else "s", 
             " not defined in model: ", paste(extra, collapse=", "))
      miss = setdiff(ng, fit)
      if (length(miss) > 0)
        stop("missing initial guess", if(length(miss)==1) "" else "es", ": ", 
             paste(miss, collapse=", "))
      private$optimizer(guess, formula, fixed, ...)
    }
  )
)
