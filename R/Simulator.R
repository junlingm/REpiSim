# ==============================================================================
# Simulation infrastructure
# ==============================================================================
#
# This file defines:
#   - default.functions: a whitelist of base-R functions that are always allowed
#     in model expressions without explicit attachment.
#   - attached.functions: an environment holding additional user-supplied
#     functions that models may call (via attach.function()).
#   - Simulator: an abstract R6 superclass for numerical simulators.
#
# Key design point:
#   Model expressions are *symbolic*. When a simulator is constructed, we verify
#   that every function referenced by the model is either:
#     (a) in default.functions, or
#     (b) present in attached.functions
#   This helps catch misspelled or missing helper functions early.
#
# IMPORTANT:
#   attached.functions MUST inherit from baseenv() so that base functions
#   (e.g., exp, log, sin) remain visible when evaluating generated code.
#
# NEW (2025-12):
#   Some simulators (e.g., Gillespie backends) return only state variables.
#   Users may still request substitutions (derived variables) in `vars`.
#   We therefore compute any requested substitutions that are missing from the
#   simulator output, using the model's substitution expressions and the
#   returned trajectories.
# ==============================================================================

# ------------------------------------------------------------------------------
# Default function whitelist
# ------------------------------------------------------------------------------

## Default functions that are always available
default.functions <- c(
  "exp", "log", "log10", "sqrt", "sin", "cos", "tan", "asin", "acos", "atan",
  "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
  "abs", "sign", "round", "floor", "ceiling", "trunc",
  "min", "max", "sum", "prod", "mean", "var", "sd",
  "cumsum", "cumprod", "diff", "lag", "lead",
  "which", "is.na", "is.finite", "is.infinite", "is.nan",
  "ifelse", "all", "any", "length", "nrow", "ncol",
  "matrix", "eigen", "qr", "svd", "chol", "solve", "det",
  "rep", "seq", "seq_along", "seq_len", "replicate",
  "sample", "sort", "order", "rank",
  "rexp", "rnorm", "runif", "rpois", "rbinom", "rt", "rchisq", "rf",
  "rweibull", "rgamma", "rnbinom", "rgeom", "rhyper", "rlogis",
  "pexp", "pnorm", "punif", "ppois", "pbinom", "pt", "pchisq", "pf",
  "pweibull", "pgamma", "pnbinom", "pgeom", "phyper", "plogis",
  "qexp", "qnorm", "qunif", "qpois", "qbinom", "qt", "qchisq", "qf",
  "qweibull", "qgamma", "qnbinom", "qgeom", "qhyper", "qlogis"
)

# ------------------------------------------------------------------------------
# Attached function registry
# ------------------------------------------------------------------------------

## Functions that can be used in models in addition to default.functions.
##
## We set parent = baseenv() so base functions remain visible during evaluation
## of generated model code, even if we set a function's environment to this env.
attached.functions <- new.env(parent = asNamespace("stats"))

#' Attach functions to be used by models
#'
#' @param ... Named arguments, each specifying a function.
#' @param functions A named list of functions (alternative input form).
#'
#' @details
#' - All entries must be named.
#' - Values must be functions; a NULL value removes the function from the
#'   registry.
#' - Attached functions are stored in an environment so that generated model
#'   code can evaluate them by name.
#'
#' @export
attach.function <- function(..., functions = NULL) {
  fs <- c(list(...), functions)
  ns <- names(fs)
  if (is.null(ns) || "" %in% ns)
    stop("functions must be named")

  check <- vapply(fs, function(f) is.function(f) || is.null(f), logical(1))
  if (!all(check))
    stop(paste(ns[!check], collapse = ", "), " must be functions (or NULL to remove)")

  for (n in ns) {
    attached.functions[[n]] <- fs[[n]]  # assigning NULL removes entry
  }
}

# ------------------------------------------------------------------------------
# Simulator base class (abstract)
# ------------------------------------------------------------------------------

#' An R6 class that is the superclass for numerical simulators
#'
#' Subclasses must implement the private methods:
#' - build(model): prepare simulation code and store it in private$.model
#' - .simulate(t, y0, parms, ...): run the numerical simulation
#'
#' The main public method is `simulate()`.
#'
#' @docType class
#' @export
Simulator <- R6Class(
  "Simulator",

  private = list(
    # a representation of the model (subclass-specific)
    .model = NULL,

    # optional compiled backend artifacts (subclass-specific)
    .compiled = NULL,

    # model compartments (character vector)
    compartments = NULL,

    # indexed compartment dimensions used to flatten structured y0
    compartment_dimensions = NULL,

    # model parameters (character vector)
    parameters = NULL,

    # indexed parameter dimensions used to validate positional strata indexing
    parameter_dimensions = NULL,

    # model substitutions represented as assignment calls (dependency-sorted)
    alias = NULL,

    # ------------------------------------------------------------------------
    # Abstract hooks for subclasses
    # ------------------------------------------------------------------------

    build = function(model) NULL,

    compile = function(model, r_model) {
      stop(
        "compiled backend is not available for ",
        class(self)[[1]],
        call. = FALSE
      )
    },

    .simulate = function(t, y0, parms, ...) NULL,

    # ------------------------------------------------------------------------
    # Helpers for generating assignment code (used by some subclasses)
    # ------------------------------------------------------------------------

    format.equation = function(eq) {
      var <- eq[[2]]
      if (is.call(var)) {
        if (var[[1]] != "'")
          stop("Invalid equation ", eq)
        var <- as.name(paste0(".d.", var[[2]]))
      }
      call("<-", var, eq[[3]])
    },

    format.var = function(S, name) {
      lapply(S, function(var) {
        call("<-", as.name(var), call("[[", as.name(name), match(var, S)))
      })
    },

    format.substitution = function() {
      c(
        private$format.var(private$compartments, "y"),
        private$format.var(private$parameters, "parms"),
        lapply(private$alias, private$format.equation)
      )
    },

    # ------------------------------------------------------------------------
    # Evaluate substitutions on a trajectory table
    # ------------------------------------------------------------------------

    # Compute requested substitutions row-wise and append them to `data`.
    compute_substitutions = function(data, parms, want) {
      if (length(want) == 0) return(data)

      defined <- intersect(want, names(private$alias))
      if (length(defined) == 0) return(data)

      # Environment that can see attached functions.
      env <- new.env(parent = attached.functions)

      # Bind parameters once (constants across rows)
      if (!is.null(parms) && length(private$parameters) > 0) {
        for (p in private$parameters) assign(p, parms[[p]], envir = env)
      }

      # Evaluate row-by-row (robust; optimize later if needed)
      for (i in seq_len(nrow(data))) {
        # bind compartments
        for (c in private$compartments) {
          if (!c %in% names(data))
            stop("simulation output is missing compartment column: ", c)
          assign(c, data[[c]][[i]], envir = env)
        }

        # evaluate substitutions in dependency order; write back the requested ones
        for (nm in names(private$alias)) {
          expr <- private$alias[[nm]][[3]]
          assign(nm, eval(expr, envir = env), envir = env)

          if (nm %in% defined) {
            if (!nm %in% names(data)) data[[nm]] <- NA_real_
            data[[nm]][[i]] <- get(nm, envir = env)
          }
        }
      }

      data
    },

    validate_parameter_strata = function(parms) {
      if (length(private$parameter_dimensions) == 0) return()

      for (parameter in names(private$parameter_dimensions)) {
        dims <- private$parameter_dimensions[[parameter]]
        value <- parms[[parameter]]
        if (is.null(value)) {
          flat_names <- strata_flat_dimension_names(parameter, dims)
          if (all(flat_names %in% names(parms))) next
        }
        expected_rank <- length(dims)

        if (expected_rank == 1) {
          if (is.null(dim(value)) && !(length(value) %in% c(1, length(dims[[1]]))))
            stop("parameter ", parameter, " length must be 1 or match model strata")
          actual <- if (is.null(dim(value))) names(value) else dimnames(value)[[1]]
          private$validate_strata_names(parameter, 1, actual, dims[[1]], require_names = FALSE)
        } else {
          actual_dim <- dim(value)
          if (is.null(actual_dim) || length(actual_dim) < expected_rank)
            stop("parameter ", parameter, " must have ", expected_rank, " dimensions")
          expected_dim <- vapply(dims, length, integer(1))
          if (!identical(as.integer(actual_dim[seq_len(expected_rank)]), expected_dim))
            stop("parameter ", parameter, " dimensions must match model strata")

          actual_dimnames <- dimnames(value)
          if (!is.null(actual_dimnames) && length(actual_dimnames) >= expected_rank) {
            for (i in seq_len(expected_rank)) {
              private$validate_strata_names(parameter, i, actual_dimnames[[i]], dims[[i]], require_names = FALSE)
            }
          }
        }
      }
    },

    flatten_parameters = function(parms) {
      if (!is.list(parms) || is.data.frame(parms)) return(parms)

      np <- names(parms)
      if (is.null(np) || any(np == ""))
        stop("parms must be named")

      out <- list()
      for (name in np) {
        dims <- private$parameter_dimensions[[name]]
        value <- parms[[name]]

        if (is.null(dims)) {
          if (length(value) != 1)
            stop("parameter ", name, " length must be 1")
          out[[name]] <- value
        } else {
          flat_names <- strata_flat_dimension_names(name, dims)
          if (length(dims) == 1 && is.null(dim(value))) {
            values <- if (length(value) == 1) rep(value, length(dims[[1]])) else as.list(value)
            names(values) <- flat_names
            out <- c(out, values)
          } else {
            values <- as.list(as.vector(value))
            names(values) <- flat_names
            out <- c(out, values)
          }
        }
      }

      unlist(out)
    },

    validate_strata_names = function(parameter, dimension, actual, expected, require_names = TRUE) {
      if (is.null(actual)) {
        if (!require_names) return()
        stop(
          "parameter ", parameter, " dimension ", dimension,
          " must be named for positional stratified indexing"
        )
      }
      if (!identical(as.character(actual), as.character(expected))) {
        stop(
          "parameter ", parameter, " dimension ", dimension,
          " names must be in model strata order: ",
          paste(expected, collapse = ", ")
        )
      }
    },

    flatten_initial_values = function(y0) {
      if (!is.list(y0) || is.data.frame(y0)) return(y0)

      ny <- names(y0)
      if (is.null(ny) || any(ny == ""))
        stop("structured y0 must be a named list")

      out <- numeric()
      for (name in ny) {
        dims <- private$compartment_dimensions[[name]]
        value <- y0[[name]]

        if (is.null(dims)) {
          if (length(value) != 1)
            stop("initial value ", name, " is not a stratified compartment")
          out[[name]] <- value
        } else {
          out <- c(out, private$flatten_stratified_initial_value(name, value, dims))
        }
      }
      out
    },

    flatten_stratified_initial_value = function(name, value, dims) {
      expected_rank <- length(dims)
      if (expected_rank == 1) {
        private$validate_strata_names(name, 1, names(value), dims[[1]])
        out <- as.numeric(value[dims[[1]]])
        names(out) <- vapply(dims[[1]], function(v) strata_flat_name(name, v), character(1))
        return(out)
      }

      actual_dim <- dim(value)
      if (is.null(actual_dim) || length(actual_dim) < expected_rank)
        stop("initial value ", name, " must have ", expected_rank, " dimensions")

      expected_dim <- vapply(dims, length, integer(1))
      if (!identical(as.integer(actual_dim[seq_len(expected_rank)]), expected_dim))
        stop("initial value ", name, " dimensions must match model strata")

      actual_dimnames <- dimnames(value)
      if (is.null(actual_dimnames) || length(actual_dimnames) < expected_rank)
        stop("initial value ", name, " must have dimnames")
      for (i in seq_len(expected_rank))
        private$validate_strata_names(name, i, actual_dimnames[[i]], dims[[i]])

      grid <- expand.grid(dims, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      out <- vapply(seq_len(nrow(grid)), function(row) {
        key <- as.character(unlist(grid[row, , drop = TRUE], use.names = FALSE))
        do.call("[", c(list(value), as.list(key)))
      }, numeric(1))
      names(out) <- vapply(seq_len(nrow(grid)), function(row) {
        key <- as.character(unlist(grid[row, , drop = TRUE], use.names = FALSE))
        strata_flat_name(name, key)
      }, character(1))
      out
    }
  ),

  public = list(
    initialize = function(model, compile = FALSE) {
      fs <- model$functions
      if (!is.null(fs) && length(fs) > 0) {
        ok <- vapply(
          fs,
          function(f) f %in% default.functions || !is.null(attached.functions[[f]]),
          logical(1)
        )
        if (!all(ok)) {
          missing <- fs[!ok]
          stop("missing functions: ", paste(missing, collapse = ", "))
        }
      }

      private$compartments <- model$flat_compartments()
      private$compartment_dimensions <- model$compartment_dimensions()
      private$parameters <- model$flat_parameters()
      private$parameter_dimensions <- model$parameter_dimensions()

      subst <- model$substitutions
      private$alias <- list()
      for (n in names(subst)) {
        s <- subst[[n]]
        private$alias[[n]] <- list(as.name("<-"), as.name(n), s)
        if (!is.null(attr(s, "compartment")))
          attr(private$alias[[n]], "compartment") <- TRUE
      }

      private$.model <- private$build(model)
      try({
        environment(private$.model) <- attached.functions
      }, silent = TRUE)

      private$.compiled <- if (isTRUE(compile)) {
        private$compile(model, private$.model)
      } else {
        NULL
      }
    },

    simulate = function(t, y0, parms = NULL, vars = names(y0), ...) {
      y0 <- private$flatten_initial_values(y0)
      if (is.null(vars)) vars <- names(y0)

      ny <- names(y0)
      if (is.null(ny) || any(ny == ""))
        stop("y0 must be named")

      y0 <- y0[private$compartments]
      missing <- which(is.na(y0))
      if (length(missing) > 0)
        stop("missing initial values for ", paste(private$compartments[missing], collapse = ", "))

      extra <- setdiff(ny, private$compartments)
      if (length(extra) > 0)
        stop("extra initial values for ", paste(extra, collapse = ", "))

      if (length(private$parameters) == 0) {
        parms <- NULL
      } else {
        if (is.null(parms))
          stop("parms must be provided for parameters: ", paste(private$parameters, collapse = ", "))

        np <- names(parms)
        if (is.null(np) || any(np == ""))
          stop("parms must be named")

        private$validate_parameter_strata(parms)
        parms <- private$flatten_parameters(parms)
        np <- names(parms)

        parms <- parms[private$parameters]
        missing <- which(if (is.list(parms)) vapply(parms, is.null, logical(1)) else is.na(parms))
        if (length(missing) > 0)
          stop("missing parameter values for ", paste(private$parameters[missing], collapse = ", "))

        extra <- setdiff(np, private$parameters)
        if (length(extra) > 0)
          stop("extra parameter values for ", paste(extra, collapse = ", "))
      }

      data <- private$.simulate(t, y0, parms, ...)

      if (!is.null(colnames(data)) && colnames(data)[[1]] == "time")
        vars <- unique(c("time", vars))

      # Compute substitutions if they were requested but not returned by backend.
      missing_vars <- setdiff(vars, names(data))
      if (length(missing_vars) > 0) {
        data <- private$compute_substitutions(data, parms, missing_vars)
      }

      still_missing <- setdiff(vars, names(data))
      if (length(still_missing) > 0) {
        stop(
          "requested vars not found in output: ",
          paste(still_missing, collapse = ", "),
          ". Available columns are: ",
          paste(names(data), collapse = ", ")
        )
      }

      data[, vars, drop = FALSE]
    }
  ),

  active = list(
    model = function() private$.model,

    compiled = function() private$.compiled
  )
)
