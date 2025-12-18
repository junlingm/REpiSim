# ==============================================================================
# Bayesian calibration base class
# ==============================================================================
#
# This file defines the Bayesian calibrator superclass.
#
# Bayesian is a thin wrapper around Calibrator that:
#   - chooses a default simulator backend (ODE) suitable for likelihood-based fitting
#   - enforces that `guess` and `priors` refer to the same set of fitted quantities
#   - forwards priors to subclass implementations via Calibrator's `...`
#
# Subclasses (e.g., Metrop) implement the actual sampling algorithm in
# private$.calibrate().
# ==============================================================================

#' A calibrator that uses a Bayesian method
#'
#' @name Bayesian
#' @docType class
#' @export
Bayesian <- R6::R6Class(
  "Bayesian",
  inherit = Calibrator,
  
  private = list(
    # Default simulator for Bayesian calibration is numerical ODE integration.
    # Subclasses can override by redefining private$simulator(model).
    simulator = function(model) {
      ODE$new(model)
    }
  ),
  
  public = list(
    #' @description
    #' Calibrate the model to data with Bayesian priors.
    #'
    #' @param initial.values Named initial values for compartments. Unspecified values
    #'   (or values set to NA / Expression, depending on your convention) are fitted.
    #' @param parms Named parameter values. Unspecified values are fitted.
    #' @param priors Named list of Distribution objects specifying priors for *fitted*
    #'   quantities. Names must match `guess`.
    #' @param guess Named numeric vector of initial guesses for all fitted quantities.
    #' @param ... Extra arguments forwarded to subclass calibrators (e.g., MCMC control).
    calibrate = function(initial.values, parms, priors, guess, ...) {
      # guess must be named
      ng <- names(guess)
      if (is.null(ng) || any(ng == ""))
        stop("guess must be a named numeric vector")
      
      # priors must be a named list of Distribution objects
      if (!is.list(priors))
        stop("priors must be a named list of Distribution objects")
      
      ok <- vapply(priors, function(x) inherits(x, "Distribution"), logical(1))
      if (!all(ok))
        stop("invalid priors for: ", paste(names(priors)[!ok], collapse = ", "))
      
      np <- names(priors)
      if (is.null(np) || any(np == ""))
        stop("priors must be named")
      
      # guess and priors must cover exactly the same fitted quantities
      if (length(setdiff(ng, np)) > 0 || length(setdiff(np, ng)) > 0) {
        stop("guess and priors must have the same names")
      }
      
      # Forward priors through ... so subclasses can receive them as an argument
      # in private$.calibrate(guess, formula, fixed, priors, ...).
      super$calibrate(initial.values, parms, guess, priors = priors, ...)
    }
  ),
  
  active = list(
    #' @field samples
    #' Samples from the posterior (subclasses should override).
    samples = function() NULL
  )
)
