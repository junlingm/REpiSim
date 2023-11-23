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
        if (!fn %in% private$.functions)
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