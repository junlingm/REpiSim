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

`%+%.Expression` <- function(a, b) {
  x = a$clone()
  x$add(b)
}

`%-%.Expression` <- function(a, b) {
  x = a$clone()
  x$sub(b)
}

`%*%.Expression` <- function(a, b) {
  x = a$clone()
  x$mul(b)
}

`%/%.Expression` <- function(a, b) {
  x = a$clone()
  x$div(b)
}

## the builtin functions that are not reported in Model$functions
builtin.functions = c(
  "+", "-", "*", "/", "^", "==", ">", "<", ">=", "<=", "!=",
  "&&", "||", "&", "|", "!", "[", "[[", "ifelse"
)

# An R6 class to extract names and functions from an expression
Expression <- R6::R6Class(
  "Expression",
  
  private = list(
    .expr = NULL,
    .parms = NULL,
    .functions = NULL,
    
    # this function extracts all parameters (names) and functions used in the 
    # given expression. 
    extract = function(expr) {
      if (is.name(expr)) {
        var = as.character(expr)
        if (!var %in% private$.parms)
          private$.parms = c(private$.parms, var)
      } else if (is.call(expr)) {
        fn = as.character(expr[[1]])
        if (!(fn %in% private$.functions) && !(fn %in% builtin.functions))
          private$.functions = c(fn, private$.functions)
        for (a in as.list(expr)[-1])
          private$extract(a)
      }
    },
    
    do.substitute = function(expr, subs) {
      if (is.name(expr)) {
        name = as.character(expr)
        new = subs[[name]]
        if (is.null(new)) expr else new
      } else if (is.call(expr)) {
        name = as.character(expr[[1]])
        new = subs[[name]]
        if (!is.null(new)) {
          if (is.character(new)) new = as.name(new)
          if (!is.name(new))
            stop("cannot replace the function name ", name, " by ", new)
          expr[[1]] = new
        }
        n = length(as.list(expr))
        for (i in 2:n)
          expr[[i]] = private$do.substitute(expr[[i]], subs)
        expr
      } else expr
    }
  ),
  
  public = list(
    #' constructor
    #' 
    #' @param expr the expression to parse
    initialize = function(expr) {
      private$.expr = expr
      private$extract(expr)
    },
    
    #' Perform a substitute
    #' 
    #' @param subs a named list of substitutions
    #' @return an expression with the substitution
    #' @details For each name in the list, the variable with that name in the
    #' expression is replaced by the corresponding value in the list.
    substitute = function(subs) {
      private$do.substitute(private$.expr, subs)
    },
    
    #' rename a parameter or a function to a new name
    #' @param from the name to rename
    #' @param to the new name
    #' @return an invisible self
    rename = function(from, to) {
      if (!is.character(from) || !is.name(from))
        stop("from must be a character or a name")
      if (!is.character(to) || !is.name(to))
        stop("to must be a character or a name")
      subs = list()
      subs[[from]] = to
      if (from %in% private$.parms) {
        private$.expr = private$do.substitute(private$.expr, subs)
        private$.parms[private$.parms == from] = to
      } else if (from %in% private$.functions) {
        private$.expr = private$do.substitute(private$.expr, subs)
        private$.functions[private$.functions == from] = to
      }
      invisible(self)
    },

    #' assign the value of an expression
    #' @param e the expression to assign
    #' @return self
    assign = function(e) {
      if (!is(e, "Expression")) e = Expression(e)
      private$.expr = e$expr
      private$.parms = e$parms
      private$.functions = e$functions
      self
    },
    
    #' add b to this expression
    #' @param b a number, a anme or a call.
    #' @return self
    add = function(b) {
      if (!is(b, "Expression")) b = Expression$new(b)
      if (is.numeric(b$expr)) {
        if (b$expr == 0) return(self)
        if (is.numeric(private$.expr)) {
          private$.expr = private$.expr + b$expr
          return(self)
        }
      }
      if (is.numeric(private$.expr) && private$.expr == 0) return (self$assign(b))
      private$.expr = as.call(list(as.name("+"), private$.expr, b$expr))
      private$.parms = union(private$.parms, b$parms)
      private$.functions = union(private$.functions, b$functions)
      self
    },
    
    
    #' subtract an expression from this one
    #' @param b the value to subtract
    #' @return self
    sub = function(b) {
      if (!is(b, "Expression")) b = Expression$new(b)
      if (is.numeric(b$expr)) {
        if (b$expr == 0) return (self)
        if (is.numeric(private$.expr)) {
          private$.expr = private$.expr - b$expr
          return(self)
        }
      }
      if (is.numeric(private$.expr) && private$.expr == 0) {
        nb = b$clone()
        return (self$assign(nb$negate()))
      }
      private$.expr = as.call(list(as.name("-"), private$.expr, b$expr))
      private$.parms = union(private$.parms, b$parms)
      private$.functions = union(private$.functions, b$functions)
      self
    },
    
    #' negate the expression
    #' @return self
    negate = function() {
      if (is.numeric(private$.expr)) {
        private$.expr = -private$.expr
      } else {
        private$.expr = if (is.call(private$.expr) && private$.expr[[1]] == "-" && length(private$.expr) == 2) {
          private$.expr[[2]]
        } else as.call(list(as.name("-"), private$.expr))
      }
      self
    },
    
    #' multiple by a factor
    #' @param b the value to multiply by
    #' @return self
    mul = function(b) {
      if (!is(b, "Expression")) b = Expression$new(b)
      if (is.numeric(b$expr)) {
        if (b$expr == 1) return(self)
        if (b$expr == 0) {
          private$.expr = 0
          private$.parms = character(0)
          private$.functions = character(0)
          return(self)
        }
        if (is.numeric(private$.expr)) {
          private$.expr = private$.expr * b$expr
          return(self)
        }
      }
      if (is.numeric(private$.expr)) {
        if (private$.expr == 0) return (self)
        if (private$.expr == 1) return (self$assign(b))
      }
      private$.expr = as.call(list(as.name("*"), private$.expr, b$expr))
      private$.parms = union(private$.parms, b$parms)
      private$.functions = union(private$.functions, b$functions)
      self
    },
    
    #' divide by a factor
    #' @param b the value to divide by
    #' @return self
    div = function(b) {
      if (is.numeric(private$.expr) && private$.expr == 0) return(self)
      if (!is(b, "Expression")) b = Expression$new(b)
      if (is.numeric(b$expr)) {
        if (b$expr == 0) stop("dividing by 0")
        if (b$expr == 1) return(self)
        if (is.numeric(private$.expr)) {
          private$.expr = private$.expr / b$expr
          return(self)
        }
      }
      private$.expr = as.call(list(as.name("/"), private$.expr, b$expr))
      private$.parms = union(private$.parms, b$parms)
      private$.functions = union(private$.functions, b$functions)
      self
    }
  ),
  
  active = list(
    #' @field parms the parameters used int he expression
    parms = function() {
      private$.parms
    },
    
    #' @field functions the function named used int he expression
    functions = function() {
      private$.functions
    },
    
    #' @field expr the original expression
    expr = function() {
      private$.expr
    }
  )
)


# build up the order of dependence of an alias
build.order = function(order, info) {
  # if var is already available in order, no need to change the order.
  if (info$name %in% order) return(order)
  # calculate the unavailable dependencies (i.e., not in order)
  deps = setdiff(private$.formula[[info$value]]$depend, order)
  # no dependencies or dependencies are already available
  # then we can just put it in order
  if (length(deps) == 0) return(c(order, info$name))
  for (d in deps) {
    v = private$.where[[d]]
    if (!is.null(v)) order = private$build.order(order, v)
  }
  c(order, info$name)
}