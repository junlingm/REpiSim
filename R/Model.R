D = function(x, ...) {
  # Helper to create a quoted expression with additional arguments.
  #
  # Example:
  #   D(S, t)  ->  `'`(S, t)
  #
  # NOTE: This helper is used elsewhere in the package; it is kept here
  # unchanged for backward compatibility.
  vars = as.list(substitute(list(...)))[-1]
  as.call(c(list("'", x), vars))
}

# ==============================================================================
# Model (R6)
# ==============================================================================
#
# The Model class represents a symbolic ODE model:
#
#   dX/dt = f(X, parameters, t)
#
# A Model is defined using R formulas:
#   S ~ -beta*S*I/N,
#.  I ~ beta*S*I/N - gamma*I,
#.  R ~ gamma*I,
#.  N = S + I + R
# representing 
#.  dS/dt = -beta*S*I/N, 
#.  dI/dt = beta*S*I/N - gamma*I,
#.  dR/dt = gamma*I,
# where N = S + I + R.
# Here:
#   - "compartments" are the state variables (e.g., S, I, R),
#   - "parameters" are the remaining free symbols appearing in a formula (e.g., beta and gamma)
#.  - In addition,  "where" substitutions can be defined for algebraic relationships (e.g., N = S + I + R),
#
# The This class is intentionally lightweight and purely symbolic: it stores
# expressions, performs name checking/renaming, tracks parameters, and exposes
# a standard representation used by downstream classes (e.g., simulators and 
# calibrators).
#
# ==============================================================================

#' R6 class representing a symbolic ODE model
#'
#' A `Model` is described by a system of ordinary differential equations (ODEs).
#' State variables are referred to as **compartments**. Each compartment is
#' defined by a formula `name ~ rate`, representing the derivative
#' \eqn{d(name)/dt = rate}.
#'
#' The model may also include algebraic substitutions (defined via `where()`),
#' which are useful for derived quantities such as a total population size
#' `N = S + I + R`. Symbols used in equations that are neither compartments nor
#' substitutions are treated as **parameters**, and are detected automatically.
#'
#' A `Model` can be:
#' - simulated numerically via the `ODE` class (wrapper around `deSolve`),
#' - simulated stochastically (for compartmental models) via `RGillespie`,
#' - typeset into LaTeX via `TexFormatter`.
#'
#' @section Construction:
#' `Model$new(...)` accepts a mixture of:
#' - formulas of the form `name ~ expr` (passed to `compartment()`), and
#' - named expressions of the form `name = expr` (passed to `where()`).
#'
#' @section Active bindings:
#' - `compartments`: character vector of compartment names.
#' - `equations`: list with two elements: `equations` (ODEs) and `where`
#'   (substitutions), both as expressions.
#' - `parameters`: character vector of parameter names.
#' - `substitutions`: named list of substitution expressions (topologically
#'   ordered by dependence; see details on the field itself).
#' - `representation`: serialisable list representation of the model.
#' - `functions`: functions referenced in any expression.
#' - `t`: name of the independent variable.
#'
#' @docType class
#' @export
Model <- R6Class(
  "Model",
  
  private = list(
    # --------------------------------------------------------------------------
    # Internal storage
    # --------------------------------------------------------------------------
    .compartments = list(),  # named list: compartment name -> list(name, value)
    .where = list(),         # named list: substitution name -> list(name, value)
    .formula = list(),       # named list: internal formula name -> Expression object
    .t = NULL,               # independent variable name (string). Invariant: ALWAYS a character scalar.
    
    # --------------------------------------------------------------------------
    # Internal helpers
    # --------------------------------------------------------------------------
    
    # Define (or redefine) an internal formula with validity checking.
    #
    # - If `formula` is NULL, the formula is removed.
    # - If a formula with the same name exists, it is overwritten.
    # - Validity checks prevent redefining compartments/substitutions/parameters
    #   as functions.
    define.formula = function(name, formula) {
      # If formula is NULL, remove the stored formula and return.
      if (is.null(formula)) {
        private$.formula[[name]] = NULL
        return()
      }
      
      e = if (is(formula, "Expression")) formula else Expression$new(formula)
      
      # Validate: do not allow redefining known names as functions.
      for (f in e$functions) {
        if (!is.null(private$.compartments[[f]]))
          stop("Redefining the compartment ", f, " as a function")
        if (!is.null(private$.where[[f]]))
          stop("Redefining the substitution ", f, " as a function")
        
        # Check for parameter collisions: if `f` is a parameter of any stored
        # formula, it cannot be used as a function name.
        for (ff in private$.formula) {
          if (f %in% ff$parms)
            stop("Redefining the parameter ", ff, " as a function")
        }
      }
      
      private$.formula[[name]] = e
    },
    
    # Perform a rename after sanity checks have been done by public$rename().
    do.rename = function(from, to) {
      if (!is.null(private$.compartments[[from]])) {
        # Rename a compartment: its derivative formula is stored under ".d.<name>".
        formula.from = paste0(".d.", from)
        formula.to = paste0(".d.", to)
        
        private$.formula[[formula.to]] = private$.formula[[formula.from]]
        private$.formula[[formula.from]] = NULL
        
        private$.compartments[[to]] = private$.compartments[[from]]
        private$.compartments[[to]]$name = to
        private$.compartments[[to]]$value = formula.to
        private$.compartments[[from]] = NULL
        
      } else if (!is.null(private$.where[[from]])) {
        # Rename a substitution: its formula is stored under its own name.
        private$.where[[to]] = private$.where[[from]]
        private$.where[[to]]$name = to
        private$.where[[to]]$value = to
        private$.where[[from]] = NULL
        
        private$.formula[[to]] = private$.formula[[from]]
        private$.formula[[from]] = NULL
        
      } else {
        stop(from, " is not defined")
      }
      
      # Apply the rename inside every stored expression.
      for (e in private$.formula) {
        e$rename(from, to)
      }
    },
    
    # Load the model from a file on disk.
    # The file is expected to load an object named `model` that is a
    # `representation` list produced by the active binding `representation`.
    load = function(file) {
      if (is.character(file)) {
        if (!file.exists(file))
          stop("file does not exist: ", file)
        e = new.env()
        load(file, envir = e)
        file = e$model
        if (is.null(file))
          stop("not a valid model file: ", file)
      }
      private$construct(file)
    },
    
    # Reconstruct the model from a representation list.
    construct = function(representation) {
      if (!identical(representation$class, "Model"))
        stop("invalid model file")
      
      private$.t = representation$.t
      # the default independent variable name is "t"
      if (!is.character(private$.t) || private$.t == "") private$.t = "t"
      
      for (C in names(representation$compartments)) {
        if (C == ".t") {
          private$.t = as.name(representation$compartments[[C]])
        } else {
          r = representation$compartments[[C]]
          self$compartment(call("~", as.name(C), r))
        }
      }
      self$where(pairs = representation$substitutions)
    }
  ),
  
  public = list(
    # --------------------------------------------------------------------------
    # Public API
    # --------------------------------------------------------------------------
    
    #' @description
    #' Construct a Model object with compartments and substitutions.
    #'
    #' @param ... A mixture of:
    #' - formulas of the form `name ~ expr` (passed to `compartment()`), and
    #' - named expressions `name = expr` (passed to `where()`).
    #' @param t Name of the independent variable, either a string (recommended)
    #' or `NULL`/empty to default to `"t"`.
    #' @param file If not `NULL`, a path (or representation object) from which
    #' to load the model.
    #'
    #' @examples
    #' # An SIR model
    #' SIR <- Model$new(
    #'   S ~ -beta*S*I/N,
    #'   I ~  beta*S*I/N - gamma*I,
    #'   R ~  gamma*I,
    #'   N = S + I + R
    #' )
    #' print(SIR)
    initialize = function(..., t = "t", file = NULL) {
      if (!is.null(file)) private$load(file)
      
      private$.t = if (is.null(t) || t == "") "t" else if (is.character(t)) t else
        stop("Invalid independent variable name ", t)
      
      args = as.list(substitute(list(...)))[-1]
      ns = names(args)
      
      if (length(args) > 0) {
        for (i in seq_along(args)) {
          if (!is.null(ns) && ns[i] != "") {
            self$where(pairs = args[i])
          } else {
            self$compartment(args[[i]])
          }
        }
      }
    },
    
    #' @description
    #' Define (or redefine) a compartment using a formula.
    #'
    #' @param eq A formula of the form `name ~ rate`. This defines a compartment
    #' with the given name and a rate of change (derivative).
    #' @return The invisible object `self` for chained operations.
    #'
    #' @details
    #' The compartment name cannot conflict with the independent variable name
    #' or an existing substitution.
    #'
    #' If `rate` is `NULL`, the compartment is deleted (equivalent to
    #' `delete(name)`).
    #'
    #' @examples
    #' # An SIR model built incrementally
    #' SIR <- Model$new()
    #' SIR$compartment(S ~ -beta*S*I/N)$
    #'   compartment(I ~  beta*S*I/N - gamma*I)$
    #'   compartment(R ~  gamma*I)$
    #'   where(N = S + I + R)
    #' print(SIR)
    compartment = function(eq) {
      if (!is.call(eq) || eq[[1]] != "~")
        stop("Invalid equation")
      
      name = as.character(eq[[2]])
      
      if (name == private$.t)
        stop(name, " is the independent variable, and so cannot be used as a dependent variable name")
      if (!is.null(private$.where[[name]]))
        stop(name, " is already defined as a substitution, so cannot be used as a dependent variable name")
      
      formula = eq[[3]]
      
      if (is.null(formula)) {
        self$delete(name)
      } else {
        rate = paste0(".d.", name)
        private$define.formula(rate, formula)
        private$.compartments[[name]] = list(name = name, value = rate)
      }
      
      invisible(self)
    },
    
    #' @description
    #' Define algebraic substitutions (named expressions).
    #'
    #' @param ... Named expressions `name = expr`.
    #' @param pairs A named list of expressions; an alternative way to provide
    #' substitutions.
    #'
    #' @details
    #' Substitutions are stored as expressions. When a substitution depends on a
    #' compartment (directly or indirectly), the returned expression in
    #' `substitutions` is marked with attribute `"compartment" = TRUE`.
    #'
    #' @return If called with no arguments, returns `self$substitutions`.
    #'
    #' @examples
    #' SIR <- Model$new()
    #' SIR$compartment(S ~ -beta*S*I)$
    #'   compartment(I ~  beta*S*I - gamma*I)$
    #'   compartment(R ~  gamma*I)$
    #'   where(pairs = list(beta = quote(b/N)))$
    #'   where(N = S + I + R)
    #' print(SIR)
    where = function(..., pairs = NULL) {
      # Combine ... and pairs (both are supported for convenience).
      defs = c(as.list(substitute(list(...)))[-1], pairs)
      
      # If called without arguments, return all substitutions.
      if (length(defs) == 0)
        return(self$substitutions)
      
      # Sanity checks
      ns = names(defs)
      if (is.null(ns) || any(ns == ""))
        stop("Substitutions must be in the form name=value")
      
      u = unique(ns)
      if (length(u) != length(ns))
        stop(
          "Redefined substitution",
          if (length(u) - length(ns) > 1) "s " else " ",
          paste(unique(setdiff(ns, u)), collapse = ", ")
        )
      
      compartments = self$compartments
      redef = ns[ns %in% compartments]
      if (length(redef) > 0)
        stop(
          "Redefining compartment",
          if (length(redef) > 1) "s " else " ",
          paste(redef, collapse = ", ")
        )
      
      # Define each substitution.
      for (name in ns) {
        if (name == private$.t)
          stop(name, " is the independent variable, and so cannot be used as a parameter name")
        private$define.formula(name, defs[[name]])
        private$.where[[name]] = list(name = name, value = name)
      }
      
      invisible(self)
    },
    
    #' @description
    #' Delete a compartment or a substitution.
    #'
    #' @param name Character (or name) specifying the compartment or substitution
    #' to be deleted. A character vector is accepted and deleted element-wise.
    #' @return The invisible object `self` for chained operations.
    #'
    #' @details
    #' Deleting a compartment removes its differential equation. Deleting a
    #' substitution removes its algebraic definition. In both cases, the name
    #' becomes an ordinary symbol; if it is still used elsewhere in the model it
    #' will appear as a parameter.
    #'
    #' @examples
    #' SIR <- Model$new(
    #'   S ~ -beta*S*I/N,
    #'   I ~  beta*S*I/N - gamma*I,
    #'   R ~  gamma*I,
    #'   N = S + I + R
    #' )
    #' SIR$delete("N")
    #' print(SIR)
    delete = function(name) {
      if (length(name) > 1) {
        for (n in name) self$delete(n)
        return(invisible(self))
      }
      
      if (!is.name(name) && !is.character(name))
        stop("invalid name ", name)
      
      if (name == private$.t)
        stop(name, " is the independent variable, and so cannot be deleted")
      
      if (!is.null(private$.compartments[[name]])) {
        private$.compartments[[name]] = NULL
        private$define.formula(paste0(".d.", name), NULL)
      } else if (!is.null(private$.where[[name]])) {
        private$.where[[name]] = NULL
        private$define.formula(name, NULL)
      } else {
        stop(name, " is not defined")
      }
      
      invisible(self)
    },
    
    #' @description
    #' Rename a compartment, a substitution, or the independent variable.
    #'
    #' @param formula A rename specification in the form `to <- from`.
    #' @param from The name to be changed (character or name).
    #' @param to The new name (character or name). If `NULL`, the target is deleted.
    #' @return The invisible object `self` for chained operations.
    #'
    #' @details
    #' The new name cannot conflict with an existing compartment, substitution,
    #' or the independent variable name.
    #'
    #' @examples
    #' SIR <- Model$new(
    #'   S ~ -beta*S*I/N,
    #'   I ~  beta*S*I/N - gamma*I,
    #'   R ~  gamma*I,
    #'   N = S + I + R
    #' )
    #' SIR$rename(U <- S)$rename(b <- beta)
    #' print(SIR)
    rename = function(formula = NULL, from = NULL, to = NULL) {
      formula = substitute(formula)
      if (!is.null(formula)) {
        if (!is.call(formula) || formula[[1]] != "<-")
          stop("invalid formula")
        to = formula[[2]]
        from = formula[[3]]
      }
      
      if (is.null(from))
        stop("The name to be changed is missing")
      if (!is.name(from) && !is.character(from))
        stop("Invalid name ", from)
      
      from = as.character(from)
      
      # Special case: renaming the independent variable
      if (from == private$.t) {
        # Renaming the independent variable.
        #
        # Invariant: private$.t is always stored as a character string, so
        # we coerce the new name (`to`) to character.
        if (is.null(to))
          stop("The independent variable cannot be deleted")
        private$.t <- as.character(to)
        return(invisible(self))
      }
      
      if (is.null(to))
        return(self$delete(from))
      if (!is.name(to) && !is.character(to))
        stop("Invalid name ", to)
      
      to = as.character(to)
      
      if (to == private$.t)
        stop(from, " cannot be renamed to the independent variable")
      
      if (from == to) return(invisible(self))
      
      private$do.rename(from, to)
      invisible(self)
    },
    
    #' @description
    #' Format the model for printing.
    format = function() {
      l = c("Model:", paste("independent variable:", private$.t))
      
      if (length(private$.compartments) > 0) {
        l = c(l, "  Compartments:")
        for (c in private$.compartments) {
          f = private$.formula[[c$value]]
          l = c(l, paste0("    ", c$name, " ~ ", deparse(f$expr)))
        }
      }
      
      if (length(private$.where) > 0) {
        l = c(l, "  where")
        for (w in private$.where) {
          l = c(l, paste0("    ", w$name, " = ", deparse(private$.formula[[w$value]]$expr)))
        }
      }
      
      if (length(self$parameters) > 0)
        l = c(l, paste0("  Parameters: ", paste(self$parameters, collapse = ", ")))
      
      paste(l, collapse = "\n")
    }
  ),
  
  active = list(
    #' @field compartments
    #' A read-only field returning a character vector of compartment names.
    compartments = function() {
      names(private$.compartments)
    },
    
    #' @field equations
    #' A read-only field returning a list with two elements:
    #' - `equations`: named list of ODE equations as expressions
    #' - `where`: named list of substitutions as expressions
    equations = function() {
      compartments = lapply(
        private$.compartments,
        function(C) {
          call("==", call("'", as.name(C$name)), private$.formula[[C$value]]$expr)
        }
      )
      where = lapply(
        private$.where,
        function(w) {
          call("==", as.name(w$name), private$.formula[[w$value]]$expr)
        }
      )
      list(equations = compartments, where = where)
    },
    
    #' @field parameters
    #' A read-only field returning a character vector of inferred parameter names.
    parameters = function() {
      parms = c()
      for (e in private$.formula) {
        if (length(e$parms) > 0) {
          parms = union(parms, e$parms)
        }
      }
      remove = c(names(private$.compartments), names(private$.where), private$.t)
      setdiff(parms, remove)
    },
    
    #' @field substitutions
    #' A read-only field returning a named list of substitution expressions.
    #'
    #' The substitutions are sorted by dependence: if substitution `A` depends on
    #' `B`, then `B` appears earlier in the list than `A`. If a substitution
    #' depends (directly or indirectly) on a compartment, the expression returned
    #' will have attribute `"compartment" = TRUE`.
    substitutions = function() {
      l = lapply(
        private$.where,
        function(w) { private$.formula[[w$value]] }
      )
      
      # The Expression class defines an ordering (dependency order).
      l = l[order(l)]
      
      nl = names(l)
      compartments = self$compartments
      
      is.compartment = sapply(
        l,
        function(e) { any(e$parms %in% compartments) }
      )
      
      # A substitution is indirectly compartment-dependent if it depends on
      # another substitution that is compartment-dependent.
      indirect = sapply(
        nl,
        function(name) {
          e = l[[name]]
          check = sapply(
            e$parms,
            function(p) { p %in% nl && is.compartment[[p]] }
          )
          any(check)
        }
      )
      
      mapply(function(e, c, i) {
        x = e$expr
        if (c || i) attr(x, "compartment") = TRUE
        x
      }, l, is.compartment, indirect)
    },
    
    #' @field representation
    #' A read-only field returning a serialisable representation of the model.
    representation = function() {
      list(
        class = "Model",
        .t = private$.t,
        compartments = sapply(
          private$.compartments,
          function(C) private$.formula[[C$value]]$expr
        ),
        substitutions = self$substitutions
      )
    },
    
    #' @field functions
    #' A read-only field returning function names referenced in expressions.
    functions = function() {
      f = c()
      for (e in private$.formula) {
        f = union(f, e$functions)
      }
      f
    },
    
    #' @field t
    #' The independent variable name (stored internally as a character string).
    t = function() {
      private$.t
    }
  )
)

if (exists("TEST") && isTRUE(TEST)) {
  m = Model$new(
    S ~ - beta*S*I/N,
    I ~ beta*S*I/N - gamma*I,
    R ~ gamma*I,
    N = S+I+R
  )
  print(m)
}
