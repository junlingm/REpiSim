# ==============================================================================
# Calibration infrastructure
# ==============================================================================
#
# This file defines the abstract Calibrator class.
#
# A calibrator fits a model to a dataset. Conceptually it needs:
#   - a model (Model / Compartmental)
#   - a simulator backend (chosen by subclasses)
#   - a dataset with observation columns
#   - a mapping from observation columns -> model variables or expressions
#   - a split of unknowns into:
#       * fixed values (provided)
#       * formula-defined values (computed from other values)
#       * fitted values (estimated by the calibrator)
#
# Subclasses must implement:
#   - private$simulator(model): construct a Simulator instance
#   - private$.calibrate(guess, formula, fixed, ...): perform optimization/inference
#   - private$interpret(results): convert algorithm output into user-facing form
# ==============================================================================

#' The base R6 class for calibrators
#'
#' @details
#' A calibrator fits a model to a dataset. Like a simulator, calibration requires:
#' - initial conditions (some fixed, some unknown),
#' - parameter values (some fixed, some unknown),
#' - data to fit,
#' - a mapping from model variables to observation variables.
#'
#' Subclasses implement the actual fitting algorithm in `private$.calibrate()`.
#'
#' @name Calibrator
#' @docType class
#' @export
Calibrator <- R6::R6Class(
  "Calibrator",
  
  private = list(
    # Simulator used during calibration (subclass-specific, e.g., ODE)
    .simulator = NULL,
    
    # Data used for fitting (numeric vector or data.frame, depending on mapping)
    .data = NULL,
    
    # A cloned model used for calibration (may include extra observation aliases)
    .model = NULL,
    
    # Whether the observed quantity is cumulative (difference is taken)
    .cumulative = FALSE,
    
    # Named character vector: observation column -> model variable name
    .mapping = NULL,
    
    # Simulation time vector corresponding to the data
    .time = NULL,
    
    # Raw output from the underlying calibration algorithm
    .details = NULL,
    
    # ------------------------------------------------------------------------
    # Hooks for subclasses
    # ------------------------------------------------------------------------
    
    #' Construct a simulator for calibration (must be implemented by subclasses)
    simulator = function(model) NULL,
    
    #' Interpret the algorithm output into a user-facing result
    interpret = function(results) NULL,
    
    #' The actual calibration routine (must be implemented by subclasses)
    .calibrate = function(guess, formula, fixed, ...) NULL,
    
    # ------------------------------------------------------------------------
    # Core helpers
    # ------------------------------------------------------------------------
    
    # Simulate the model given:
    #   pars   - named numeric vector of fitted values (subset of unknowns)
    #   formula- named list of Expression objects (computed values)
    #   fixed  - named list of fixed values
    #
    # Returns a matrix/data.frame aligned with the data being fitted.
    simulate = function(pars, formula, fixed, ...) {
      # Combine fitted and fixed into one mutable "value table".
      all <- c(as.list(pars), fixed)
      
      # Evaluate formula-defined values (dependency order assumed).
      for (n in names(formula)) {
        f <- formula[[n]]
        all[[n]] <- eval(f$expr, envir = all)
      }
      
      # Extract initial conditions and parameters
      ic <- unlist(all[private$.model$compartments])
      par <- all[private$.model$parameters]
      
      # Simulate only the variables needed by the mapping (plus time)
      vars <- unlist(private$.mapping)
      sim <- private$.simulator$simulate(private$.time, ic, par, vars = vars, ...)
      
      # If cumulative observations were provided, convert cumulative model outputs
      # to per-interval increments by differencing.
      if (private$.cumulative) {
        if (ncol(sim) == 2) {
          # single series: time + one variable
          return(diff(sim[[2]]))
        }
        return(as.data.frame(lapply(sim[, -1, drop = FALSE], diff)))
      }
      
      # Drop time column for fitting
      sim[, -1, drop = FALSE]
    },
    
    # Split a named vector/list into: fixed numeric values, formula-defined values,
    # and variables-to-fit (those not specified).
    #
    # For `mode="initial.value"`, the set is model compartments.
    # For `mode="parameter"`, the set is model parameters.
    split = function(x, mode = c("initial.value", "parameter")) {
      mode <- match.arg(mode)
      
      set <- if (mode == "initial.value") {
        private$.model$compartments
      } else {
        private$.model$parameters
      }
      
      if (length(x) == 0) {
        return(list(value = NULL, formula = NULL, fit = set))
      }
      
      ns <- names(x)
      if (is.null(ns) || any(ns == ""))
        stop(mode, " must be named")
      
      extra <- setdiff(ns, set)
      if (length(extra) > 0)
        stop(paste(extra, collapse = ", "), " not defined in the model")
      
      idx.values <- vapply(x, is.numeric, logical(1))
      idx.formula <- vapply(x, function(y) is(y, "Expression"), logical(1))
      
      x.values <- as.list(x[idx.values])
      x.formula <- x[idx.formula]
      x.fit <- setdiff(set, ns[idx.values | idx.formula])
      
      list(value = x.values, formula = x.formula, fit = x.fit)
    },
    
    # Infer what to fit, what to compute from formula, and what is fixed.
    #
    # Returns a list suitable for do.call(private$.calibrate, ...):
    #   list(guess=..., formula=..., fixed=..., <extra args>)
    fit.info = function(initial.values, parms, guess, ...) {
      ic <- private$split(initial.values, mode = "initial.value")
      p <- private$split(parms, mode = "parameter")
      
      formula <- c(ic$formula, p$formula)
      
      # Formula objects are Expression instances; the package provides an `order()`
      # helper that orders by dependence.
      if (length(formula) > 0) formula <- formula[order(formula)]
      
      # Collect additional symbols referenced by formulas that are neither
      # compartments nor parameters: these become "extra" fitted values.
      needed <- Reduce(function(extra, f) union(f$parms, extra), formula, init = character())
      extra <- setdiff(needed, c(private$.model$compartments, private$.model$parameters))
      
      fit <- c(extra, ic$fit, p$fit)
      if (length(fit) == 0)
        stop("no initial values or parameters to fit")
      
      if (!is.numeric(guess) || any(is.na(guess)))
        stop("invalid initial guess values")
      
      ng <- names(guess)
      if (is.null(ng) || any(ng == ""))
        stop("guess must be a named numeric vector")
      
      extra_guess <- setdiff(ng, fit)
      if (length(extra_guess) == 1) {
        stop("guess contains extra parameter: ", extra_guess)
      } else if (length(extra_guess) > 1) {
        stop("guess contains extra parameters: ", paste(extra_guess, collapse = ", "))
      }
      
      missed <- setdiff(fit, ng)
      if (length(missed) == 1) {
        stop("guess is missing parameter: ", missed)
      } else if (length(missed) > 1) {
        stop("guess is missing parameters: ", paste(missed, collapse = ", "))
      }
      
      fixed <- c(ic$value, p$value)
      
      c(
        list(
          guess = guess,
          formula = formula,
          fixed = fixed
        ),
        list(...)
      )
    }
  ),
  
  public = list(
    #' @description
    #' Initialize a calibrator.
    #'
    #' @param model The model to calibrate.
    #' @param time Either:
    #'   - numeric vector of times (including initial time if `cumulative=TRUE`), or
    #'   - character column name in `data` that contains time.
    #' @param data A data.frame with observation columns.
    #' @param ... Optional mapping formulas of the form `data_col ~ model_var_or_expr`.
    #' @param cumulative Whether the data are cumulative (if TRUE, the model output
    #'   is differenced before fitting).
    #' @param mapping Optional named character vector mapping data columns to model
    #'   variables (alternative to `...` formulas).
    #'
    #' @details
    #' A mapping entry maps a data column name to either:
    #' - a model compartment/substitution name, or
    #' - an expression in model variables (in which case a new substitution is
    #'   created internally).
    initialize = function(model, time, data, ..., cumulative = FALSE, mapping = character()) {
      m <- model$clone(deep = TRUE)
      
      if (!is.data.frame(data))
        stop("data must be a data.frame object")
      
      if (length(mapping) > 0) {
        extra_cols <- setdiff(names(mapping), colnames(data))
        if (length(extra_cols) > 0)
          stop(
            "data column", if (length(extra_cols) == 1) "" else "s",
            " does not exist: ", paste(extra_cols, collapse = ", ")
          )
      }
      
      private$.mapping <- mapping
      private$.cumulative <- cumulative
      
      # Resolve time vector
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
      
      # Parse mapping formulas in ...
      args <- as.list(substitute(list(...)))[-1]
      if (length(args) == 0 && length(mapping) == 0)
        stop("no mapping from data to model variables")
      
      for (a in args) {
        if (!is.call(a) || as.character(a[[1]]) != "~")
          stop("invalid mapping: ", deparse(a))
        
        # LHS: data column
        if ((!is.name(a[[2]]) && !is.character(a[[2]])) || is.null(data[[as.character(a[[2]])]]))
          stop("not a data column: ", deparse(a[[2]]))
        col <- as.character(a[[2]])
        
        # RHS: model variable name OR expression in model variables
        if (is.name(a[[3]])) {
          var <- as.character(a[[3]])
          if (!var %in% m$compartments && is.null(m$substitutions[[var]]))
            stop(var, " is not a model variable")
        } else {
          if (!is.call(a[[3]]))
            stop(deparse(a[[3]]), " is not a valid specification")
          
          # If the data column name does not already conflict, reuse it;
          # otherwise create an internal observation alias.
          var <- if (!col %in% m$compartments && is.null(m$substitutions[[col]])) {
            col
          } else {
            paste0(".observation.", col)
          }
          
          m$where(pairs = setNames(list(a[[3]]), var))
        }
        
        if (col %in% names(private$.mapping))
          stop("redefining mapping for data column ", col)
        
        private$.mapping[col] <- var
      }
      
      private$.model <- m
      private$.simulator <- private$simulator(m)
      
      # Keep only the mapped observation columns in the fitting data
      obs_cols <- names(private$.mapping)
      fit_data <- data[, obs_cols, drop = FALSE]
      if (ncol(fit_data) == 1) fit_data <- fit_data[[1]]
      private$.data <- fit_data
    },
    
    #' Calibrate the model to data
    #'
    #' @param initial.values Named initial values (compartments). Values may be:
    #'   - numeric (fixed),
    #'   - Expression (computed), or
    #'   - omitted (fit).
    #' @param parms Named parameter values; same conventions as initial.values.
    #' @param guess Named numeric vector giving starting values for all fitted quantities.
    #' @param ... Extra arguments forwarded to the subclass calibrator.
    calibrate = function(initial.values, parms, guess, ...) {
      info <- private$fit.info(initial.values, parms, guess, ...)
      private$.details <- do.call(private$.calibrate, info)
      private$interpret(private$.details)
    }
  ),
  
  active = list(
    #' @field parameters Names of all model parameters
    parameters = function() private$.model$parameters,
    
    #' @field details Raw output from the underlying calibration algorithm
    details = function() private$.details
  )
)
