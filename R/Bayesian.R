#' A calibrator that uses a Baysian method.
#' @name Baysian
#' @docType class
#' @export
#' 

Baysian <- R6::R6Class(
  "Baysian",
  inherit = Calibrator,
  
  private = list(
    simulator = function(model) {
      ODE$new(model)
    }
  ),
  
  public = list (
    #' @description Calibrate the model to data
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
      ng = names(guess)
      if (is.null(ng) || any(ng == ""))
        stop("guess must be a named list")
      if (!is.list(priors) || !all(sapply(priors, function(x) "Distribution" %in% class(x))))
        stop("invalid priors")
      np = names(priors)
      if (is.null(np) || any(np == ""))
        stop("priors must be a named list")
      # check that ng and ns have the same names
      if (length(setdiff(ng, np)) > 0 || length(setdiff(np, ng)) > 0)
        stop("guess and priors must have the same names")
      super$calibrate(initial.values, parms, guess, priors, ...)
    }
  ),
  
  active = list(
    #' @field samples the samples from the fit, a read-only field
    samples = function() {
      NULL
    }
  )
)
