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

calibrator.is.missing <- function(x) {
  is.null(x) || (length(x) == 1 && is.atomic(x) && is.na(x))
}

calibrator.is.formula <- function(x) {
  is(x, "Expression") || is.call(x) || is.name(x)
}

calibrator.as.expression <- function(x) {
  if (is(x, "Expression")) x else Expression$new(x)
}

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
    #' @field .compile Whether default simulators should use compiled backends
    .compile = FALSE,
    #' @field .mapping Named character vector mapping data columns -> model variables
    .mapping = NULL,
    #' @field .time Numeric vector of times used for simulation
    .time = NULL,
    #' @field .details Raw output returned by the subclass algorithm
    .details = NULL,
    #' @field .compartment_dimensions Stratified compartment dimensions
    .compartment_dimensions = list(),
    #' @field .parameter_dimensions Stratified parameter dimensions
    .parameter_dimensions = list(),
    
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
      par <- private$reconstruct_parameters(all)
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
    #' @param x Named list/vector of values and/or formulas.
    #' @param mode One of "initial.value" or "parameter".
    split = function(x, mode = c("initial.value", "parameter")) {
      mode <- match.arg(mode)
      x <- private$normalize_calibration_values(x, mode)
      set <- if (mode == "initial.value") private$.model$compartments else private$calibration_parameters()
      
      if (length(x) == 0) return(list(value = NULL, formula = NULL, fit = set))
      
      ns <- names(x)
      if (is.null(ns) || any(ns == "")) stop(mode, " must be named")
      
      extra <- setdiff(ns, set)
      if (length(extra) > 0) stop(paste(extra, collapse = ", "), " not defined in the model")
      
      idx.formula <- vapply(x, calibrator.is.formula, logical(1))
      idx.values <- vapply(
        seq_along(x),
        function(i) is.numeric(x[[i]]) && !idx.formula[[i]] && !calibrator.is.missing(x[[i]]),
        logical(1)
      )
      
      x.values <- as.list(x[idx.values])
      x.formula <- lapply(x[idx.formula], calibrator.as.expression)
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
      
      guess <- private$normalize_guess(guess)
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
    },

    calibration_parameters = function() {
      params <- private$.model$parameters
      expanded <- unlist(
        lapply(names(private$.parameter_dimensions), function(parameter) {
          strata_flat_dimension_names(parameter, private$.parameter_dimensions[[parameter]])
        }),
        use.names = FALSE
      )
      c(setdiff(params, c(names(private$.parameter_dimensions), expanded)), expanded)
    },

    normalize_calibration_values = function(x, mode = c("initial.value", "parameter")) {
      mode <- match.arg(mode)
      if (length(x) == 0) return(x)

      if (!is.list(x) || is.data.frame(x)) x <- as.list(x)
      ns <- names(x)
      if (is.null(ns) || any(ns == "")) return(x)

      out <- list()
      for (name in ns) {
        value <- x[[name]]
        dims <- if (mode == "initial.value") {
          private$.compartment_dimensions[[name]]
        } else {
          private$.parameter_dimensions[[name]]
        }

        if (is.null(dims)) {
          out[[name]] <- value
        } else {
          out <- c(out, private$flatten_dimensioned_value(name, value, dims, mode))
        }
      }
      out
    },

    normalize_guess = function(guess) {
      if (length(guess) == 0) return(guess)
      if (!is.list(guess) || is.data.frame(guess)) guess <- as.list(guess)

      ns <- names(guess)
      if (is.null(ns) || any(ns == "")) return(unlist(guess))

      out <- list()
      for (name in ns) {
        value <- guess[[name]]
        if (!is.null(private$.compartment_dimensions[[name]])) {
          out <- c(out, private$flatten_dimensioned_value(
            name, value, private$.compartment_dimensions[[name]], "guess"
          ))
        } else if (!is.null(private$.parameter_dimensions[[name]])) {
          out <- c(out, private$flatten_dimensioned_value(
            name, value, private$.parameter_dimensions[[name]], "guess"
          ))
        } else {
          out[[name]] <- value
        }
      }
      out <- out[!vapply(out, calibrator.is.missing, logical(1))]
      unlist(out)
    },

    flatten_dimensioned_value = function(name, value, dims, mode) {
      if (length(dims) == 1 && is.null(dim(value))) {
        if (!(length(value) %in% c(1, length(dims[[1]]))))
          stop(mode, " ", name, " length must be 1 or match model strata")
        if (identical(mode, "guess") && !is.null(names(value)) &&
            all(as.character(names(value)) %in% as.character(dims[[1]]))) {
          values <- as.list(value)
          names(values) <- vapply(names(value), function(v) strata_flat_name(name, v), character(1))
          return(values)
        }
        if (!is.null(names(value)) && !identical(as.character(names(value)), as.character(dims[[1]])))
          stop(mode, " ", name, " names must be in model strata order: ", paste(dims[[1]], collapse = ", "))

        values <- if (length(value) == 1) rep(value, length(dims[[1]])) else as.list(value)
        names(values) <- private$flat_dimension_names(name, dims)
        return(values)
      }

      actual_dim <- dim(value)
      if (is.null(actual_dim) || length(actual_dim) < length(dims))
        stop(mode, " ", name, " must have ", length(dims), " dimensions")

      expected_dim <- vapply(dims, length, integer(1))
      if (!identical(as.integer(actual_dim[seq_along(dims)]), expected_dim))
        stop(mode, " ", name, " dimensions must match model strata")

      actual_dimnames <- dimnames(value)
      if (!is.null(actual_dimnames) && length(actual_dimnames) >= length(dims)) {
        for (i in seq_along(dims)) {
          dn <- actual_dimnames[[i]]
          if (!is.null(dn) && !identical(as.character(dn), as.character(dims[[i]])))
            stop(mode, " ", name, " dimension ", i, " names must be in model strata order: ",
                 paste(dims[[i]], collapse = ", "))
        }
      }

      flat_names <- private$flat_dimension_names(name, dims)
      values <- as.list(as.vector(value))
      names(values) <- flat_names
      values
    },

    flat_dimension_names = function(name, dims) {
      strata_flat_dimension_names(name, dims)
    },

    reconstruct_parameters = function(all) {
      lapply(private$.model$parameters, function(parameter) {
        dims <- private$.parameter_dimensions[[parameter]]
        if (is.null(dims)) return(all[[parameter]])

        flat_names <- private$flat_dimension_names(parameter, dims)
        values <- unlist(all[flat_names], use.names = FALSE)
        if (length(values) == length(flat_names) && !any(is.na(values))) {
          out <- array(values, dim = vapply(dims, length, integer(1)), dimnames = dims)
          if (length(dims) == 1) out <- as.vector(out)
          names(out) <- if (length(dims) == 1) dims[[1]] else names(out)
          return(out)
        }

        all[[parameter]]
      }) |> stats::setNames(private$.model$parameters)
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
    initialize = function(model, time, data, ..., cumulative = FALSE, mapping = character(), compile = FALSE) {
      m <- model$clone(deep = TRUE)
      if (!is.data.frame(data)) stop("data must be a data.frame object")
      private$.compile <- isTRUE(compile)
      
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
      private$.compartment_dimensions <- m$compartment_dimensions()
      private$.parameter_dimensions <- m$parameter_dimensions()
      private$.simulator <- private$simulator(m)
      
      # PATCH: always keep observation data as a data.frame
      obs_cols <- names(private$.mapping)
      private$.data <- data[, obs_cols, drop = FALSE]
    },
    
    #' @description
    #' Calibrate the model to data.
    #'
    #' @param initial.values Named initial values; omitted values are estimated.
    #'   Values may also be quoted expressions, e.g. `I0 = quote(exp(log_I0))`.
    #' @param parms Named parameter values; omitted values are estimated. Values
    #'   may also be quoted expressions, e.g. `a = quote(exp(b))`, to fit a
    #'   transformed parameter.
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
