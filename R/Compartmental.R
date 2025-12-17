# ==============================================================================
# Compartmental (R6)
# ==============================================================================
#
# The Compartmental class is a subclass of Model. Instead of defining ODEs
# directly (e.g., S ~ -beta*S*I/N), a compartmental model is defined by
# transitions of the form:
#
#   from -> to ~ rate
#
# Each transition contributes:
#   - a negative flow to the 'from' compartment
#   - a positive flow to the 'to' compartment
#
# The ODE for each compartment is the sum of all incoming rates minus the sum of
# all outgoing rates. Because transitions correspond to discrete events, the
# same transition specification can also be used for stochastic simulation
# (e.g., Gillespie algorithms).
#
# Notes on conventions used in this implementation:
#   - Compartments must be defined before they can be used in transitions.
#   - `NULL -> X` represents an external input into X (e.g., births).
#   - `X -> NULL` represents an external output from X (e.g., deaths).
#   - The independent variable name is stored in private$.t as a STRING, to
#     match the invariant in Model (see Model.R).
#   - A per-capita transition rate can be specified either by:
#       * percapita(rate_expr) inside the formula, or
#       * percapita=TRUE with a rate that is per-capita.
#     In both cases, the total rate is multiplied by the 'from' compartment size.
#
# ==============================================================================

#' An R6 class representing a compartmental model
#'
#' `Compartmental` is a subclass of [Model]. While a `Model` defines the ODE
#' system directly, a `Compartmental` model defines **transitions** between
#' compartments. The ODE system is derived automatically by summing inflows and
#' outflows implied by those transitions.
#'
#' Because transitions correspond to events, compartmental models can also be
#' used for stochastic simulation (e.g., Gillespie methods).
#'
#' @name Compartmental
#' @docType class
#' @export
Compartmental <- R6Class(
  "Compartmental",
  inherit = Model,
  
  private = list(
    # --------------------------------------------------------------------------
    # Internal storage
    # --------------------------------------------------------------------------
    
    # Named list of transitions. Each transition record is a list with fields:
    #   - from: character or NULL
    #   - to  : character or NULL
    #   - name: character (unique key)
    .transitions = list(),
    
    # --------------------------------------------------------------------------
    # Internal helpers
    # --------------------------------------------------------------------------
    
    # (Re)formulate the ODE for a compartment by aggregating transition rates.
    #
    # For each transition:
    #   - subtract its rate from the source compartment (if non-NULL)
    #   - add its rate to the destination compartment (if non-NULL)
    #
    # This method updates the underlying Model's compartment equation.
    equation = function(compartment) {
      rate = Reduce(
        function(rate, tr) {
          r = private$.formula[[tr$name]]
          if (!is.null(tr$from) && identical(tr$from, compartment)) {
            rate %-% r
          } else if (!is.null(tr$to) && identical(tr$to, compartment)) {
            rate %+% r
          } else {
            rate
          }
        },
        private$.transitions,
        Expression$new(0)
      )
      super$compartment(call("~", as.name(compartment), rate))
    },
    
    # Parse the "side" of a transition formula, allowing an optional rate.
    #
    # Expected forms:
    #   - X            (compartment only)
    #   - X ~ expr     (compartment with rate)
    #
    # Returns:
    #   - NULL if e is NULL
    #   - a name/symbol if only a compartment is provided
    #   - list(compartment=<name>, rate=<Expression>) if a rate is included
    parse.rate = function(e) {
      if (is.null(e)) return(NULL)
      if (!is.call(e)) return(e)
      
      if (e[[1]] != "~") stop("invalid transition")
      list(compartment = e[[2]], rate = Expression$new(e[[3]]))
    },
    
    # Generate a unique name for an unnamed transition.
    transition.name = function(from, to) {
      base = paste0(from, "->", to)
      name = base
      i = 1
      while (!is.null(private$.transitions[[name]])) {
        name = paste0(base, ".", i)
        i = i + 1
      }
      name
    },
    
    # Override rename to support renaming transitions as well as model symbols.
    #
    # - If `from` matches a transition name, rename that transition.
    # - Otherwise delegate to Model's rename behavior.
    do.rename = function(from, to) {
      if (!is.null(private$.transitions[[from]])) {
        private$.transitions[[to]] = private$.transitions[[from]]
        private$.transitions[[to]]$name = to
        private$.transitions[[from]] = NULL
        
        private$.formula[[to]] = private$.formula[[from]]
        private$.formula[[from]] = NULL
      } else {
        super$do.rename(from, to)
      }
    },
    
    # Reconstruct the model from a serialisable representation.
    #
    # NOTE: The original version had typos (`privte$`, `represnetation$`) that
    # would prevent reconstruction. This implementation keeps behavior but fixes
    # those misspellings.
    construct = function(representation) {
      if (!identical(representation$class, "Compartmental"))
        stop("invalid model file")
      
      # Independent variable name: keep as string, default to "t"
      private$.t <- if (is.null(representation$.t) || representation$.t == "") {
        "t"
      } else {
        as.character(representation$.t)
      }
      
      # Compartments
      if (!is.null(representation$compartments)) {
        sapply(representation$compartments, self$compartment)
      }
      
      # Transitions
      if (!is.null(representation$transitions)) {
        sapply(representation$transitions, function(tr) do.call(self$transition, tr))
      }
      
      # Substitutions
      self$where(pairs = representation$substitutions)
    }
  ),
  
  public = list(
    #' @description
    #' Construct a compartmental model.
    #'
    #' @param ... Compartment names and/or substitutions (passed to [Model] constructor).
    #' @param t Independent variable name (string recommended).
    #' @param file If not NULL, a file/representation to load the model from.
    #'
    #' @examples
    #' SIR <- Compartmental$new(S, I, R)
    #' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
    #' SIR$transition(I->R ~ gamma*I, name="recovery")
    #' print(SIR)
    initialize = function(..., t = "t", file = NULL) {
      super$initialize(..., t = t, file = file)
    },
    
    #' @description
    #' Define a compartment (state variable) by name.
    #'
    #' This differs from [Model]$compartment(), which defines an ODE formula.
    #' Here we create a compartment with a default derivative 0, to be later
    #' filled in by transitions.
    #'
    #' @param name The compartment name (character or name).
    #' @return The invisible object `self` for chained operations.
    compartment = function(name) {
      if (length(name) == 0) return(invisible(self))
      if (!is.name(name) && !is.character(name))
        stop("invalid compartment name: ", name)
      
      name_chr = as.character(name)
      if (name_chr == private$.t)
        stop(name_chr, " is the independent variable, so cannot be used as a compartment name")
      
      super$compartment(call("~", as.name(name_chr), 0))
      invisible(self)
    },
    
    #' @description
    #' Delete a compartment, substitution, or transition.
    #'
    #' If a compartment is deleted, any transitions incident to that compartment
    #' are deleted first.
    #'
    #' @param name Name of a compartment, substitution, or transition.
    #' @return The invisible object `self` for chained operations.
    delete = function(name) {
      if (!is.null(private$.transitions[[name]])) {
        # Delete a transition by name.
        self$transition(name = name)
      } else {
        # If deleting a compartment, also remove incident transitions.
        if (!is.null(private$.compartments[[name]])) {
          for (tr in private$.transitions) {
            if (identical(tr$from, name) || identical(tr$to, name))
              self$transition(name = tr$name)
          }
        }
        super$delete(name)
      }
      invisible(self)
    },
    
    #' @description
    #' Define, modify, or delete a transition.
    #'
    #' @param formula A transition specified as `from -> to ~ rate` (entered in
    #' practice as `to <- from ~ rate` due to R syntax).
    #' @param ... Named substitutions local to this call (forwarded to `where()`).
    #' @param from Source compartment (name/character) or NULL.
    #' @param to Destination compartment (name/character) or NULL.
    #' @param rate Transition rate expression. If NULL, the transition is deleted.
    #' @param percapita If TRUE, interpret `rate` as per-capita and multiply by `from`.
    #' @param name Transition name. If NULL, one is generated.
    #'
    #' @return The transition name (invisibly). If deleted, returns NULL (invisibly).
    #'
    #' @details
    #' - Inputs: `NULL -> X` (entered as `X <- NULL ~ rate` when using `formula`).
    #' - Outputs: `X -> NULL` (entered as `NULL <- X ~ rate` when using `formula`).
    #' - Per-capita rate can be specified as `percapita(expr)` or `percapita=TRUE`.
    transition = function(formula = NULL, ..., from = NULL, to = NULL, rate = NULL, percapita = FALSE, name = NULL) {
      formula = substitute(formula)
      
      # ----------------------------------------------------------------------
      # Parse transition specification
      # ----------------------------------------------------------------------
      if (!is.null(formula)) {
        if (!is.call(formula) || formula[[1]] != "<-")
          stop("invalid transition")
        
        # Left side is "to" (may include a rate via "~")
        to_parsed = private$parse.rate(formula[[2]])
        if (is.list(to_parsed)) {
          rate_obj = to_parsed$rate
          to = to_parsed$compartment
        } else {
          rate_obj = NULL
          to = to_parsed
        }
        
        # Right side is "from" (may include rate via "~" if not on left)
        from_parsed = private$parse.rate(formula[[3]])
        if (is.list(from_parsed)) {
          if (!is.null(rate_obj)) stop("invalid transition: rate specified twice")
          rate_obj = from_parsed$rate
          from = from_parsed$compartment
        } else {
          from = from_parsed
        }
        
        if (!is.null(rate_obj)) rate = rate_obj
      } else {
        # If no formula, interpret `rate` argument as an expression unless NULL.
        if (!is.null(rate) && !is(rate, "Expression"))
          rate = Expression$new(substitute(rate))
      }
      
      # ----------------------------------------------------------------------
      # Validate endpoints and normalize to character (or NULL)
      # ----------------------------------------------------------------------
      if (!is.null(from) && !is.character(from)) {
        if (!is.name(from)) stop("invalid compartment ", from)
        from = as.character(from)
      }
      if (!is.null(to) && !is.character(to)) {
        if (!is.name(to)) stop("invalid compartment ", to)
        to = as.character(to)
      }
      
      # independent variable name cannot be used as compartment
      if (!is.null(from) && identical(from, private$.t))
        stop(from, " is the independent variable, and cannot be used as a compartment")
      if (!is.null(to) && identical(to, private$.t))
        stop(to, " is the independent variable, and cannot be used as a compartment")
      
      # compartments must exist if non-NULL
      if (!is.null(from) && is.null(private$.compartments[[from]]))
        stop("the compartment ", from, " is not defined")
      if (!is.null(to) && is.null(private$.compartments[[to]]))
        stop("the compartment ", to, " is not defined")
      
      # ----------------------------------------------------------------------
      # Interpret per-capita wrappers / flags
      # ----------------------------------------------------------------------
      if (is(rate, "Expression") && is.call(rate$expr) && identical(rate$expr[[1]], as.name("percapita"))) {
        rate = Expression$new(rate$expr[[2]])
        percapita = TRUE
      }
      
      if (isTRUE(percapita)) {
        if (is.null(from)) stop("percapita=TRUE requires a non-NULL 'from' compartment")
        # Total rate = per-capita rate * population in 'from'
        rate$mul(as.name(from))
      }
      
      # ----------------------------------------------------------------------
      # Name handling
      # ----------------------------------------------------------------------
      if (is.null(name)) name = private$transition.name(from, to)
      
      # ----------------------------------------------------------------------
      # Create / update / delete transition
      # ----------------------------------------------------------------------
      tr = private$.transitions[[name]]
      
      if (!is.null(tr)) {
        # Existing transition: either delete or modify
        if (is.null(rate) || (is.null(from) && is.null(to))) {
          # Delete
          private$define.formula(name, NULL)
          private$.transitions[[name]] = NULL
          name = NULL
        } else {
          # Modify: decide whether to auto-rename if name is based on "from->to"
          rename.auto = grepl(paste0(from, "->", to), name) &&
            (!identical(tr$from, from) || !identical(tr$to, to))
          
          new.name = if (rename.auto) private$transition.name(from, to) else name
          
          # Remove old formula and define new
          private$define.formula(name, NULL)
          private$define.formula(new.name, rate)
          
          # Update record (use new.name)
          private$.transitions[[name]] = NULL
          private$.transitions[[new.name]] = list(from = from, to = to, name = new.name)
          
          # Recompute equations for affected compartments (old and new endpoints)
          if (!is.null(tr$from)) private$equation(tr$from)
          if (!is.null(tr$to)) private$equation(tr$to)
          if (!is.null(from)) private$equation(from)
          if (!is.null(to)) private$equation(to)
          
          name = new.name
        }
      } else {
        # New transition
        if (is.null(rate) || (is.null(from) && is.null(to))) {
          # nothing to add
          name = NULL
        } else {
          private$define.formula(name, rate)
          private$.transitions[[name]] = list(from = from, to = to, name = name)
          if (!is.null(from)) private$equation(from)
          if (!is.null(to)) private$equation(to)
        }
      }
      
      # ----------------------------------------------------------------------
      # Local substitutions passed via ...
      # ----------------------------------------------------------------------
      where = as.list(substitute(list(...)))[-1]
      if (length(where) > 0) {
        ns = names(where)
        for (i in seq_along(where)) {
          w = where[i]
          n = ns[i]
          if (n == private$.t)
            stop(n, " is the independent variable, and thus cannot be redefined")
          if (is.null(n) || n == "")
            stop("meaningless definition ", w)
          if (!is.null(private$.where[[n]]))
            stop("redefinition of ", n)
          self$where(pairs = w)
        }
      }
      
      invisible(name)
    },
    
    #' @description
    #' Format the class for printing.
    format = function() {
      l = "Compartmental:"
      
      if (length(private$.compartments) > 0) {
        l = c(l, paste0("  Compartments: ", paste(self$compartments, collapse = ", ")))
      }
      
      if (length(private$.transitions) > 0) {
        l = c(l, "  transitions:")
        for (tr in private$.transitions) {
          l = c(l, paste0(
            "    \"", tr$name, "\" : ",
            if (is.null(tr$from)) "NULL" else tr$from, " -> ",
            if (is.null(tr$to)) "NULL" else tr$to,
            " ~ ", deparse(private$.formula[[tr$name]]$expr)
          ))
        }
      }
      
      if (length(private$.where) > 0) {
        l = c(l, "  where")
        s = self$substitutions
        for (w in names(s))
          l = c(l, paste0("    ", w, " = ", deparse(s[[w]])))
      }
      
      if (length(self$parameters) > 0)
        l = c(l, paste0("  Parameters: ", paste(self$parameters, collapse = ", ")))
      
      paste(l, collapse = "\n")
    }
  ),
  
  active = list(
    #' @field transitions
    #' A read-only field returning the transitions, including the evaluated rate expression.
    transitions = function() {
      lapply(
        private$.transitions,
        function(x) {
          x$rate = private$.formula[[x$name]]$expr
          x
        }
      )
    },
    
    #' @field representation
    #' A read-only field returning a serialisable representation of the model.
    #'
    #' The representation contains compartments, transitions, substitutions, and `.t`.
    representation = function() {
      list(
        class = "Compartmental",
        .t = private$.t,
        compartments = self$compartments,
        transitions = self$transitions,
        substitutions = self$substitutions
      )
    }
  )
)

if (exists("TEST") && is.logical(TEST) && TEST) {
  m = Compartmental$new(S, I, R)
  m$transition(I <- S ~ beta*S*I/N, N = S + I + R)
  m$transition(R <- I ~ gamma, percapita = TRUE)
  print(m)
}
