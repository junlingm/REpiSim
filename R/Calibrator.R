# ==============================================================================
# Calibration infrastructure
# ==============================================================================
#
# This file defines the Calibrator base class.
#
# PATCH NOTE (2025-12-17):
#   private$.data is now ALWAYS stored as a data.frame (even when only a single
#   observation column is mapped). This avoids downstream dimension errors such as:
#     "private$.data[, i] invalid number of dimensions"
#   in Bayesian/MCMC objectives that iterate over columns.
# ==============================================================================

#' Base R6 class for calibrators
#'
#' @details
#' A calibrator fits a model to a dataset. It combines:
#' - a *simulator* (e.g., ODE) to produce model solutions,
#' - a *mapping* from data columns to model variables (or expressions of them),
#' - a specification of which initial conditions / parameters are fixed, derived
#'   from formulas, or estimated.
#'
#' The base class implements:
#' - parsing/validation of the mapping
#' - splitting initial values and parameters into fixed/formula/fitted sets
#' - a generic [calibrate()] entry point that delegates to subclass algorithms
#'
#' Subclasses should implement:
#' - `private$simulator(model)` to choose a simulation backend
#' - `private$.calibrate(guess, formula, fixed, ...)` to run the fitting algorithm
#' - `private$interpret(results)` to convert raw results to a user-facing object
#'
#' @name Calibrator
#' @docType class
#' @export
Calibrator <- R6::R6Class(
  "Calibrator",
  private = list(
    #' @field .simulator Internal simulator instance used during calibration
    .simulator = NULL,
    #' @field .data Observed data (always stored as a data.frame)
    .data = NULL,
    #' @field .model Model instance used during calibration (may include observation substitutions)
    .model = NULL,
    #' @field .cumulative Whether observed data are cumulative counts/curves
    .cumulative = FALSE,
    #' @field .mapping Named character vector mapping data columns -> model variables
    .mapping = NULL,
    #' @field .time Numeric vector of times used for simulation
    .time = NULL,
    #' @field .details Raw output returned by the subclass algorithm
    .details = NULL,
    
    # ------------------------------------------------------------------------
    # Hooks for subclasses
    # ------------------------------------------------------------------------
    
    #' Create a simulator for the (possibly augmented) model.
    #' Subclasses MUST override this.
    simulator = function(model) NULL,
    
    # ------------------------------------------------------------------------
    # Core simulation helper used by many calibrators
    # ------------------------------------------------------------------------
    #
    # Given a combined list `pars` containing fitted quantities (as list entries),
    # plus:
    # - `formula`: named list of Expression objects to evaluate derived quantities
    # - `fixed`: named list of fixed values
    #
    # This method:
    # 1) evaluates all formulas in dependency order,
    # 2) builds compartment ICs and parameter values,
    # 3) runs the simulator,
    # 4) returns model outputs aligned with the observation mapping.
    simulate = function(pars, formula, fixed, ...) {
      all <- c(as.list(pars), fixed)
      for (n in names(formula)) {
        f <- formula[[n]]
        all[[n]] <- eval(f$expr, envir = all)
      }
      
      ic <- unlist(all[private$.model$compartments])
      par <- all[private$.model$parameters]
      vars <- unlist(private$.mapping)
      
      sim <- private$.simulator$simulate(private$.time, ic, par, vars = vars, ...)
      
      if (private$.cumulative) {
        # Convert cumulative solutions to incident counts via first differences.
        if (ncol(sim) == 2) return(diff(sim[[2]]))
        return(as.data.frame(lapply(sim[, -1, drop = FALSE], diff)))
      }
      
      sim[, -1, drop = FALSE]
    },
    
    #' Interpret raw calibration results.
    #' Subclasses SHOULD override this.
    interpret = function(results) NULL,
    
    #' The algorithm-specific calibration routine.
    #' Subclasses MUST override this.
    .calibrate = function(guess, formula, fixed, ...) NULL,
    
    # ------------------------------------------------------------------------
    # Utilities: split fixed/formula/fitted sets
    # ------------------------------------------------------------------------
    
    #' Split initial values or parameters into fixed / formula / fitted parts.
    #'
    #' @param x Named list/vector of values and/or Expression objects.
    #' @param mode One of "initial.value" or "parameter".
    split = function(x, mode = c("initial.value", "parameter")) {
      mode <- match.arg(mode)
      set <- if (mode == "initial.value") private$.model$compartments else private$.model$parameters
      
      if (length(x) == 0) return(list(value = NULL, formula = NULL, fit = set))
      
      ns <- names(x)
      if (is.null(ns) || any(ns == "")) stop(mode, " must be named")
      
      extra <- setdiff(ns, set)
      if (length(extra) > 0) stop(paste(extra, collapse = ", "), " not defined in the model")
      
      idx.values <- vapply(x, is.numeric, logical(1))
      idx.formula <- vapply(x, function(y) is(y, "Expression"), logical(1))
      
      x.values <- as.list(x[idx.values])
      x.formula <- x[idx.formula]
      x.fit <- setdiff(set, ns[idx.values | idx.formula])
      
      list(value = x.values, formula = x.formula, fit = x.fit)
    },
    
    #' Infer what needs to be fitted vs fixed vs computed by formula.
    #'
    #' Returns a list containing:
    #' - guess: numeric vector (algorithm-specific; can be a list in future extensions)
    #' - formula: named list of Expression objects
    #' - fixed: named list of fixed values
    #' - plus any extra arguments passed through ...
    fit.info = function(initial.values, parms, guess, ...) {
      ic <- private$split(initial.values, mode = "initial.value")
      p <- private$split(parms, mode = "parameter")
      
      formula <- c(ic$formula, p$formula)
      if (length(formula) > 0) formula <- formula[order(formula)]
      
      needed <- Reduce(function(extra, f) union(f$parms, extra), formula, init = character())
      extra <- setdiff(needed, c(private$.model$compartments, private$.model$parameters))
      
      fit <- c(extra, ic$fit, p$fit)
      if (length(fit) == 0) stop("no initial values or parameters to fit")
      
      if (!is.numeric(guess) || any(is.na(guess))) stop("invalid initial values")
      
      ng <- names(guess)
      if (is.null(ng) || any(ng == "")) stop("guess must be named")
      
      extra_guess <- setdiff(ng, fit)
      if (length(extra_guess) == 1) stop("guess contains extra parameter: ", extra_guess)
      if (length(extra_guess) > 1) stop("guess contains extra parameters: ", paste(extra_guess, collapse = ", "))
      
      missed <- setdiff(fit, ng)
      if (length(missed) == 1) stop("guess is missing parameter: ", missed)
      if (length(missed) > 1) stop("guess is missing parameters: ", paste(missed, collapse = ", "))
      
      fixed <- c(ic$value, p$value)
      c(list(guess = guess, formula = formula, fixed = fixed), list(...))
    }
  ),
  
  public = list(
    #' @description
    #' Initialize a calibrator.
    #'
    #' @param model A `Model` (or subclass) to calibrate.
    #' @param time Either (i) a numeric vector of simulation times, or (ii) a character
    #'   giving the name of the time column in `data`.
    #' @param data A data.frame containing observations.
    #' @param ... Mapping formulas of the form `data_column ~ model_variable_or_expression`.
    #' @param cumulative Logical; whether `data` represents cumulative curves.
    #' @param mapping Optional named character vector mapping data columns to model variables.
    #'
    #' @details
    #' Mapping can be specified by:
    #' - `...` formulas: e.g. `cases ~ I` or `cases ~ rho * I`
    #' - `mapping`: a named character vector, e.g. `c(cases = "I")`
    #'
    #' When the RHS of a mapping formula is an expression, the expression is added as a
    #' substitution to the internal model clone so it can be simulated like other variables.
    initialize = function(model, time, data, ..., cumulative = FALSE, mapping = character()) {
      m <- model$clone(deep = TRUE)
      if (!is.data.frame(data)) stop("data must be a data.frame object")
      
      if (length(mapping) > 0) {
        extra_cols <- setdiff(names(mapping), colnames(data))
        if (length(extra_cols) > 0)
          stop("data column", if (length(extra_cols) == 1) "" else "s",
               " does not exist: ", paste(extra_cols, collapse = ", "))
      }
      
      private$.mapping <- mapping
      private$.cumulative <- cumulative
      
      # Time specification
      if (is.numeric(time)) {
        if (cumulative) {
          if (length(time) - 1 != nrow(data))
            stop("for cumulative data, time must have one more point than nrow(data)")
        } else {
          if (length(time) != nrow(data))
            stop("time should have the same length as nrow(data)")
        }
        private$.time <- time
      } else if (is.character(time) && !is.null(data[[time]])) {
        if (cumulative)
          stop("for cumulative data, provide a numeric time vector with one extra point")
        private$.time <- data[[time]]
      } else {
        stop("time must be either a numeric vector or a column name in `data`")
      }
      
      # Parse mapping formulas from ...
      args <- as.list(substitute(list(...)))[-1]
      if (length(args) == 0 && length(mapping) == 0)
        stop("no mapping from data to model variables")
      
      for (a in args) {
        if (!is.call(a) || as.character(a[[1]]) != "~")
          stop("invalid mapping: ", deparse(a))
        
        if ((!is.name(a[[2]]) && !is.character(a[[2]])) || is.null(data[[as.character(a[[2]])]]))
          stop("not a data column: ", deparse(a[[2]]))
        col <- as.character(a[[2]])
        
        if (is.name(a[[3]])) {
          var <- as.character(a[[3]])
          if (!var %in% m$compartments && is.null(m$substitutions[[var]]))
            stop(var, " is not a model variable")
        } else {
          if (!is.call(a[[3]]))
            stop(deparse(a[[3]]), " is not a valid specification")
          var <- if (!col %in% m$compartments && is.null(m$substitutions[[col]])) col else paste0(".observation.", col)
          m$where(pairs = setNames(list(a[[3]]), var))
        }
        
        if (col %in% names(private$.mapping))
          stop("redefining mapping for data column ", col)
        
        private$.mapping[col] <- var
      }
      
      private$.model <- m
      private$.simulator <- private$simulator(m)
      
      # PATCH: always keep observation data as a data.frame
      obs_cols <- names(private$.mapping)
      private$.data <- data[, obs_cols, drop = FALSE]
    },
    
    #' @description
    #' Calibrate the model to data.
    #'
    #' @param initial.values Named initial values; values to estimate are indicated
    #'   by omission and/or algorithm-specific conventions.
    #' @param parms Named parameter values; values to estimate are indicated similarly.
    #' @param guess Initial guesses for fitted quantities (format may be subclass-specific).
    #' @param ... Additional arguments forwarded to the subclass algorithm.
    calibrate = function(initial.values, parms, guess, ...) {
      info <- private$fit.info(initial.values, parms, guess, ...)
      private$.details <- do.call(private$.calibrate, info)
      private$interpret(private$.details)
    }
  ),
  
  active = list(
    #' @field parameters Names of model parameters (read-only)
    parameters = function() private$.model$parameters,
    #' @field details Raw output of the underlying algorithm (read-only)
    details = function() private$.details
  )
)
