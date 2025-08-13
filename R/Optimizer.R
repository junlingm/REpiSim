#' A subclass of Calibrator that uses an optimizer and an initial guess.
#' @name Optimizer
#' @docType class
#' @export
Optimizer <- R6::R6Class(
  "Optimizer",
  inherit = Calibrator,
  
  private = list(
    simulator = function(model) {
      ODE$new(model)
    },
    
    ## The objective function for optimization
    ## @param par the values of the fitted parameters
    ## @param formula gives the parameter values that should be calculated
    ## using this formula
    ## @param fixed the parameter values that are given
    ## @return the fitting results, which will be passed to the intepret method.
    ## @details This method must be implemented by subclasses.
    objective = function(par, formula, fixed, ...) {
      stop("The objective method must be implemented by a subclass")
    },
    
    ## The optimization is done in this function.
    ## This method must be implemented by subclasses.
    ## @param guess the initial guesses of the fitted parameters
    ## @param formula gives the parameter values that should be calculated using this formula
    ## @param fixed the parameter values that are given
    ## @return the fitting results, which will be passed to 
    ## the intepret method.
    optimizer = function(guess, formula, fixed, ...) {
      stop("The optimizer method must be implemented by a subclass")
    },
    
    ## The actual calibration is done in this function.
    ## This method Must be implemented by subclasses.
    ## @param guess the initial guess of the parameters to be fitted
    ## @param formula gives the parameter values that should be calculated using this formula
    ## @param fixed the parameter values that are given
    ## @return the fitting results, which will be passed to 
    ## the intepret method.
    .calibrate = function(guess, formula, fixed, ...) {
      private$optimizer(guess, formula, fixed, ...)
    }
  )
)
