# ==============================================================================
# ODE simulator (deSolve backend)
# ==============================================================================
#
# This file defines:
#   - ODE: a Simulator subclass that compiles a Model into a function compatible
#          with deSolve::ode().
#   - Equilibrium: a helper subclass that searches for an equilibrium by
#          repeated forward simulation ("fixed-point" iteration).
#
# Key implementation idea for ODE:
#   deSolve expects a function with signature:
#     f(time, state, parameters) -> list(dstate_dt, extra_outputs)
#
# We therefore compile:
#   1) assignments that bind y and parms into named variables
#   2) assignments for substitution variables (aliases)
#   3) derivative expressions for each compartment
#   4) return list(c(dS, dI, ...), c(extra1, extra2, ...))
#
# The "extra outputs" allow deSolve to carry any compartment-dependent
# substitutions along with the solution (so that simulate() can return them
# without recomputing).
# ==============================================================================

#' A R6 class for numerically solving ODE models
#'
#' `ODE` is a subclass of [Simulator]. This class is available only if the
#' **deSolve** package is installed.
#'
#' @name ODE
#' @docType class
#' @export
ODE <- R6Class(
  "ODE",
  inherit = Simulator,
  
  private = list(
    header = "#include \"cpp11.hpp\"
#include <cmath>
",

    setup = "
cpp11::environment attached_functions(R_NilValue);

[[cpp11::register]]
void REpiSimSetup(cpp11::environment functions)
{
  attached_functions = functions;
}
",

    # --------------------------------------------------------------------------
    # Compile model -> deSolve function(time, y, parms)
    # --------------------------------------------------------------------------
    build = function(model) {
      system <- model$equations
      
      # Derivative variable names: ".d.S", ".d.I", ...
      der <- lapply(names(system$equations), function(var) as.name(paste0(".d.", var)))
      
      # Extra outputs (substitutions): return these from the deSolve function so
      # they are included in the output data.frame.
      alias_names <- names(private$alias)
      extra_syms <- lapply(alias_names, as.name)
      names(extra_syms) <- alias_names
      
      # Function body:
      # {
      #   <bind y/parms>
      #   <compute aliases>
      #   <compute derivatives>
      #   list(c(dS, dI, ...), c(extra1, extra2, ...))
      # }
      body_calls <- c(
        as.name("{"),
        private$format.substitution(),
        lapply(system$equations, private$format.equation),
        call(
          "list",
          as.call(c(list(as.name("c")), der)),
          as.call(c(list(as.name("c")), extra_syms))
        )
      )
      
      # Build function arguments: (t, y, parms) but with t named as model$t
      args <- alist(,,)
      names(args) <- c(model$t, "y", "parms")
      
      as.function(c(args, as.call(body_calls)))
    },

    build_rhs_cpp = function(model) {
      cpp <- CppConverter$new(private$compartments, private$parameters, private$alias)
      system <- model$equations

      n_der <- length(system$equations)
      n_extra <- length(private$alias)
      out <- "__out"
      der_names <- paste0("__der_", seq_len(n_der) - 1)

      deriv_stmts <- lapply(seq_along(system$equations), function(i) {
        cpp$assign(der_names[[i]], cpp$expr(system$equations[[i]][[3]]), "double")
      })

      result <- c(
        paste0("cpp11::writable::doubles ", out, "(", n_der + n_extra, ");"),
        lapply(seq_len(n_der), function(i) {
          cpp$assign(cpp$array(out, i - 1), der_names[[i]])
        }),
        lapply(seq_along(private$alias), function(i) {
          nm <- names(private$alias)[[i]]
          cpp$assign(cpp$array(out, n_der + i - 1), cpp$expr(as.name(nm)))
        }),
        paste0("return ", out, ";")
      )

      stmts <- c(
        "cpp11::doubles __y = y;",
        "cpp11::doubles __parms = parms;",
        cpp$bind_vars(private$compartments, "y"),
        cpp$bind_vars(private$parameters, "parms"),
        cpp$attached_functions(),
        cpp$substitutions(),
        deriv_stmts,
        result
      )

      paste0(
        "[[cpp11::register]]\n",
        "cpp11::doubles ode_rhs(double ", model$t,
        ", cpp11::doubles y, cpp11::doubles parms) ",
        cpp$block(stmts, "")
      )
    },

    compile = function(model, r_model) {
      if (!requireNamespace("cpp11", quietly = TRUE)) {
        stop(
          "cpp11 is required to compile ODE. Use compile = FALSE or install cpp11",
          call. = FALSE
        )
      }

      rhs <- private$build_rhs_cpp(model)
      program <- paste(private$header, private$setup, rhs, sep = "\n")

      lib <- cpp11::cpp_source(code = program)
      rhs_fun <- get("ode_rhs", inherits = TRUE)
      REpiSimSetup(attached.functions)

      n_der <- length(private$compartments)
      extra_names <- names(private$alias)
      private$.model <- private$compiled_model(model, rhs_fun, n_der, extra_names)

      list(
        program = program,
        rhs = rhs,
        lib = lib,
        evaluate = rhs_fun
      )
    },

    compiled_model = function(model, rhs_fun, n_der, extra_names) {
      args <- alist(,,)
      names(args) <- c(model$t, "y", "parms")

      f <- function() {}
      body(f) <- substitute({
        if (is.null(parms)) parms <- numeric()
        values <- rhs_fun(time, y, parms)
        dy <- values[seq_len(n_der)]
        extra <- values[-seq_len(n_der)]

        if (length(extra_names) > 0) {
          names(extra) <- extra_names
        } else {
          extra <- numeric()
        }

        list(dy, extra)
      }, list(time = as.name(model$t)))

      formals(f) <- args
      environment(f) <- environment()
      f
    },
    
    # --------------------------------------------------------------------------
    # Run deSolve
    # --------------------------------------------------------------------------
    .simulate = function(t, y0, parms, ...) {
      as.data.frame(deSolve::ode(y = y0, times = t, func = self$model, parms = parms, ...))
    },

    finalize = function() {
      lib <- private$.compiled$lib
      if (!is.null(lib) && !is.null(lib[["path"]])) {
        try(dyn.unload(lib[["path"]]), silent = TRUE)
      }
    }
  ),
  
  public = list(
    #' @description
    #' Construct an ODE simulator for a model.
    #'
    #' @param model A [Model] (or subclass) to simulate.
    initialize = function(model, compile = FALSE) {
      if (!requireNamespace("deSolve", quietly = TRUE))
        stop("package 'deSolve' is not installed; the ODE simulator is unavailable")
      super$initialize(model, compile = compile)
    }
  ),
  
  active = list(
    #' @field model
    #' A read-only field returning the compiled deSolve function used as `func`.
    model = function() {
      private$.model
    }
  )
)

# ==============================================================================
# Equilibrium helper
# ==============================================================================

#' An R6 class to approximate an equilibrium of a model
#'
#' This class searches for an equilibrium by repeated forward simulation:
#' simulate forward, take the final state, use it as the next initial condition,
#' and stop when the 2-norm difference falls below `error`.
#'
#' @param model A [Model] object.
#' @param max.iter Maximum number of iterations.
#' @param vary Optional parameter name to vary (character).
#' @param range Numeric vector of parameter values used when `vary` is not NULL.
#' @param error Convergence tolerance on the 2-norm.
#' @param vars Variables (compartments/substitutions) to return.
#'
#' @name Equilibrium
#' @docType class
#' @export
Equilibrium <- R6Class(
  "Equilibrium",
  inherit = ODE,
  
  private = list(
    vary = NULL,
    range = NULL,
    max.iter = NULL,
    error = NULL,
    
    # Compute equilibrium for one parameter setting by fixed-point iteration
    eq = function(t, y0, parms, ...) {
      y_prev <- y0
      
      for (i in seq_len(private$max.iter)) {
        # We want the full final state each iteration, so request states only.
        sol <- super$simulate(t, y_prev, parms, vars = names(y_prev), ...)
        y_last <- unlist(sol[nrow(sol), -1, drop = TRUE])
        
        e <- sqrt(sum((y_prev - y_last)^2))
        if (!is.na(e) && e < private$error)
          return(sol[nrow(sol), -1, drop = FALSE])
        
        y_prev <- y_last
      }
      
      # Non-convergence: return NA row with correct column structure
      y_prev[] <- NA_real_
      as.data.frame(as.list(y_prev))
    },
    
    .simulate = function(t, y0, parms, ...) {
      out <- data.frame()
      
      if (is.null(private$vary)) {
        out <- rbind(out, private$eq(t, y0, parms, ...))
      } else {
        for (v in private$range) {
          parms[[private$vary]] <- v
          out <- rbind(out, private$eq(t, y0, parms, ...))
        }
      }
      
      out
    }
  ),
  
  public = list(
    #' @description
    #' Construct an equilibrium finder.
    #'
    #' @param model A [Model] (or subclass).
    #' @param vary Optional parameter name to vary.
    #' @param range Numeric vector giving the values of `vary`.
    #' @param max.iter Maximum iterations for the fixed-point iteration.
    #' @param error Convergence tolerance.
    initialize = function(model, vary = NULL, range = NULL, max.iter = 100, error = 1e-6) {
      super$initialize(model)
      private$vary <- vary
      private$range <- range
      
      if (!is.null(vary)) {
        if (!is.character(vary) || !(vary %in% model$parameters))
          stop("vary must be a parameter name")
        if (!is.numeric(range))
          stop("range must be a numeric vector")
      }
      
      private$max.iter <- max.iter
      private$error <- error
    },
    
    #' @description
    #' Compute equilibria and return them as a data frame.
    #'
    #' If `vary` is set, the output includes an extra column with the parameter
    #' values in `range`.
    simulate = function(t, y0, parms = NULL, vars = names(y0), ...) {
      data <- super$simulate(t, y0, parms, vars, ...)
      
      if (!is.null(private$vary)) {
        col <- list()
        col[[private$vary]] <- private$range
        data <- cbind(as.data.frame(col), data)
      }
      
      data
    }
  )
)
