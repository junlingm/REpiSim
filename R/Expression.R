# ==============================================================================
# Expression utilities
# ==============================================================================
#
# This file provides:
#   1) Generic arithmetic operators for Expression objects: %+%  %-%
#      %*%  %/% (defined as S3 generics and dispatched on class "Expression").
#   2) The Expression R6 class, which:
#        - stores an R language object (call/name/numeric) in $expr
#        - extracts symbol names ("parameters") and non-builtin function names
#        - supports safe symbolic renaming and substitution
#        - supports basic symbolic arithmetic (add/sub/mul/div) with small
#          simplifications (e.g., x + 0 = x, x * 1 = x, x * 0 = 0).
#   3) A dependency-order helper for lists of Expression objects.
#
# NOTE on naming:
#   - The helper `order()` is intentionally kept for backward compatibility,
#     but it masks base::order if attached on the search path. Prefer calling it
#     with explicit namespace (e.g., REpiSim::order) or using `expr_order()`.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# S3 generics for Expression arithmetic
# ------------------------------------------------------------------------------

#' Define a generic addition operator
#' @param a the first argument
#' @param b the second argument
#' @return the sum of a and b
#' @export
`%+%` <- function(a, b) UseMethod("%+%", a)

#' Define a generic subtraction operator
#' @param a the first argument
#' @param b the second argument
#' @return the difference of a and b
#' @export
`%-%` <- function(a, b) UseMethod("%-%", a)

#' Define a generic multiplication operator
#' @param a the first argument
#' @param b the second argument
#' @return the product of a and b
#' @export
`%*%` <- function(a, b) UseMethod("%*%", a)

#' Define a generic division operator
#' @param a the first argument
#' @param b the second argument
#' @return the quotient of a and b
#' @export
`%/%` <- function(a, b) UseMethod("%/%", a)

#' @export
`%+%.Expression` <- function(a, b) {
  x <- a$clone(deep = TRUE)
  x$add(b)
}

#' @export
`%-%.Expression` <- function(a, b) {
  x <- a$clone(deep = TRUE)
  x$sub(b)
}

#' @export
`%*%.Expression` <- function(a, b) {
  x <- a$clone(deep = TRUE)
  x$mul(b)
}

#' @export
`%/%.Expression` <- function(a, b) {
  x <- a$clone(deep = TRUE)
  x$div(b)
}

# ------------------------------------------------------------------------------
# Internal: builtin functions
# ------------------------------------------------------------------------------

# Built-in functions/operators that should NOT be reported in Model$functions.
# This list is conservative; you can extend it if new internal operators appear.
builtin.functions <- c(
  "+", "-", "*", "/", "^", "==", ">", "<", ">=", "<=", "!=",
  "&&", "||", "&", "|", "!", "[", "[[", "ifelse", "("
)

# ------------------------------------------------------------------------------
# Expression class
# ------------------------------------------------------------------------------

#' An R6 class to extract names and functions from an expression
#'
#' The `Expression` class wraps an R language object (name/call/numeric, etc.)
#' and provides:
#' - `$parms`: character vector of symbol names appearing in the expression
#' - `$functions`: character vector of non-builtin function names called
#' - `$substitute(subs)`: syntactic substitution on symbols and function names
#' - `$rename(from, to)`: rename a symbol or function name
#' - basic symbolic arithmetic with lightweight simplification
#'
#' @docType class
#' @export
Expression <- R6::R6Class(
  "Expression",
  
  private = list(
    .expr = NULL,
    .parms = character(0),
    .functions = character(0),
    
    # Extract parameter names and non-builtin function names from expr.
    #
    # - Names (symbols) are recorded in `.parms`.
    # - Calls record their head in `.functions` unless builtin.
    extract = function(expr) {
      if (is.name(expr)) {
        var <- as.character(expr)
        if (!(var %in% private$.parms))
          private$.parms <- c(private$.parms, var)
        
      } else if (is.call(expr)) {
        fn <- as.character(expr[[1]])
        if (!(fn %in% private$.functions) && !(fn %in% builtin.functions))
          private$.functions <- c(fn, private$.functions)
        
        for (a in as.list(expr)[-1])
          private$extract(a)
      }
    },
    
    # Recursive syntactic substitution.
    #
    # - If expr is a name and subs contains that name, replace it.
    # - If expr is a call and subs contains the function name, replace the head
    #   *only if the replacement is a valid name*.
    # - Recurse into call arguments.
    do.substitute = function(expr, subs) {
      if (is.name(expr)) {
        name <- as.character(expr)
        new <- subs[[name]]
        if (is.null(new)) expr else new
        
      } else if (is.call(expr)) {
        name <- as.character(expr[[1]])
        
        # Substitute function name if requested (must remain a name).
        new <- subs[[name]]
        if (!is.null(new)) {
          if (is.character(new)) new <- as.name(new)
          if (!is.name(new))
            stop("cannot replace the function name ", name, " by ", new)
          expr[[1]] <- new
        }
        
        # Substitute call arguments.
        n <- length(expr)
        if (n >= 2) {
          for (i in 2:n)
            expr[[i]] <- private$do.substitute(expr[[i]], subs)
        }
        expr
        
      } else {
        expr
      }
    }
  ),
  
  public = list(
    #' @description
    #' Construct an Expression object and extract symbols/functions.
    #'
    #' @param expr A language object (name/call/numeric), typically produced by
    #' `quote()` or `substitute()`.
    initialize = function(expr) {
      private$.expr <- expr
      private$.parms <- character(0)
      private$.functions <- character(0)
      private$extract(expr)
    },
    
    #' Perform syntactic substitution
    #'
    #' @param subs A named list of substitutions. Values should be language
    #' objects (e.g., `as.name("x")`, `quote(a+b)`) or character (for renaming
    #' function heads).
    #' @return A language object representing the substituted expression.
    substitute = function(subs) {
      private$do.substitute(private$.expr, subs)
    },
    
    #' Rename a symbol or a function name
    #'
    #' @param from Name to rename (character or name).
    #' @param to New name (character or name).
    #' @return The invisible object `self` for chained operations.
    rename = function(from, to) {
      if (!is.character(from) && !is.name(from))
        stop("from must be a character or a name")
      if (!is.character(to) && !is.name(to))
        stop("to must be a character or a name")
      
      from_chr <- as.character(from)
      to_chr <- as.character(to)
      
      subs <- list()
      subs[[from_chr]] <- to_chr
      
      if (from_chr %in% private$.parms) {
        private$.expr <- private$do.substitute(private$.expr, subs)
        private$.parms[private$.parms == from_chr] <- to_chr
      } else if (from_chr %in% private$.functions) {
        private$.expr <- private$do.substitute(private$.expr, subs)
        private$.functions[private$.functions == from_chr] <- to_chr
      }
      
      invisible(self)
    },
    
    #' Assign from another Expression (or raw expression)
    #'
    #' @param e Another `Expression` or a raw language object.
    #' @return `self`.
    assign = function(e) {
      if (!is(e, "Expression")) e <- Expression$new(e)
      private$.expr <- e$expr
      private$.parms <- e$parms
      private$.functions <- e$functions
      self
    },
    
    #' Add to this expression (in-place)
    #'
    #' @param b A number, name, call, or `Expression`.
    #' @return `self`.
    add = function(b) {
      if (!is(b, "Expression")) b <- Expression$new(b)
      
      # Simplify numeric additions
      if (is.numeric(b$expr)) {
        if (b$expr == 0) return(self)
        if (is.numeric(private$.expr)) {
          private$.expr <- private$.expr + b$expr
          return(self)
        }
      }
      if (is.numeric(private$.expr) && private$.expr == 0)
        return(self$assign(b))
      
      private$.expr <- as.call(list(as.name("+"), private$.expr, b$expr))
      private$.parms <- union(private$.parms, b$parms)
      private$.functions <- union(private$.functions, b$functions)
      self
    },
    
    #' Subtract from this expression (in-place)
    #'
    #' @param b The value to subtract (number/name/call/Expression).
    #' @return `self`.
    sub = function(b) {
      if (!is(b, "Expression")) b <- Expression$new(b)
      
      # Simplify numeric subtraction
      if (is.numeric(b$expr)) {
        if (b$expr == 0) return(self)
        if (is.numeric(private$.expr)) {
          private$.expr <- private$.expr - b$expr
          return(self)
        }
      }
      
      # 0 - b = -(b)
      if (is.numeric(private$.expr) && private$.expr == 0) {
        nb <- b$clone(deep = TRUE)
        return(self$assign(nb$negate()))
      }
      
      private$.expr <- as.call(list(as.name("-"), private$.expr, b$expr))
      private$.parms <- union(private$.parms, b$parms)
      private$.functions <- union(private$.functions, b$functions)
      self
    },
    
    #' Negate this expression (in-place)
    #'
    #' @return `self`.
    negate = function() {
      if (is.numeric(private$.expr)) {
        private$.expr <- -private$.expr
      } else {
        # Simplify -(-x) = x
        private$.expr <- if (is.call(private$.expr) &&
                             identical(private$.expr[[1]], as.name("-")) &&
                             length(private$.expr) == 2) {
          private$.expr[[2]]
        } else {
          as.call(list(as.name("-"), private$.expr))
        }
      }
      self
    },
    
    #' Multiply this expression (in-place)
    #'
    #' @param b The value to multiply by (number/name/call/Expression).
    #' @return `self`.
    mul = function(b) {
      if (!is(b, "Expression")) b <- Expression$new(b)
      
      # Simplify numeric multiplication
      if (is.numeric(b$expr)) {
        if (b$expr == 1) return(self)
        if (b$expr == 0) {
          private$.expr <- 0
          private$.parms <- character(0)
          private$.functions <- character(0)
          return(self)
        }
        if (is.numeric(private$.expr)) {
          private$.expr <- private$.expr * b$expr
          return(self)
        }
      }
      
      if (is.numeric(private$.expr)) {
        if (private$.expr == 0) return(self)
        if (private$.expr == 1) return(self$assign(b))
      }
      
      private$.expr <- as.call(list(as.name("*"), private$.expr, b$expr))
      private$.parms <- union(private$.parms, b$parms)
      private$.functions <- union(private$.functions, b$functions)
      self
    },
    
    #' Divide this expression (in-place)
    #'
    #' @param b The value to divide by (number/name/call/Expression).
    #' @return `self`.
    div = function(b) {
      if (is.numeric(private$.expr) && private$.expr == 0) return(self)
      if (!is(b, "Expression")) b <- Expression$new(b)
      
      # Simplify numeric division
      if (is.numeric(b$expr)) {
        if (b$expr == 0) stop("dividing by 0")
        if (b$expr == 1) return(self)
        if (is.numeric(private$.expr)) {
          private$.expr <- private$.expr / b$expr
          return(self)
        }
      }
      
      private$.expr <- as.call(list(as.name("/"), private$.expr, b$expr))
      private$.parms <- union(private$.parms, b$parms)
      private$.functions <- union(private$.functions, b$functions)
      self
    }
  ),
  
  active = list(
    #' @field parms
    #' The parameter (symbol) names used in the expression.
    parms = function() {
      private$.parms
    },
    
    #' @field functions
    #' The non-builtin function names used in the expression.
    functions = function() {
      private$.functions
    },
    
    #' @field expr
    #' The underlying language object.
    expr = function() {
      private$.expr
    }
  )
)

# ------------------------------------------------------------------------------
# Dependency ordering
# ------------------------------------------------------------------------------

#' Return the dependency order of expressions
#'
#' @param exprs A named list of `Expression` objects.
#' @return A character vector giving the order in which expressions should be
#' evaluated so that dependencies are available.
#'
#' @details
#' If expression `A` uses symbol `B` and `B` is also provided in `exprs`, then
#' `B` appears earlier than `A` in the returned order.
#'
#' NOTE: This function name masks `base::order` when attached; consider using
#' [expr_order()] instead for clarity.
#' @export
order <- function(exprs) {
  if (length(exprs) == 0) return(NULL)
  
  ns <- names(exprs)
  if (is.null(ns) || any(ns == ""))
    stop("expressions must be named")
  
  # The exprs must be a list of Expression objects
  check <- vapply(exprs, function(e) is(e, "Expression"), logical(1))
  if (!all(check))
    stop(paste(ns[!check], collapse = ", "), " must be Expression objects")
  
  build <- function(out, name) {
    if (name %in% out) return(out)
    e <- exprs[[name]]
    for (p in e$parms) {
      if (!is.null(exprs[[p]]))
        out <- build(out, p)
    }
    c(out, name)
  }
  
  Reduce(build, ns, init = character(0))
}

#' Return the dependency order of expressions (alias)
#'
#' @param exprs A named list of `Expression` objects.
#' @return Same as [order()].
#' @export
expr_order <- function(exprs) {
  order(exprs)
}
