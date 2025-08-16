## the default functions that are always available
default.functions = c(
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
  "rweibull", "rgamma", "rnbinom", "rgeom", "rhyper", "rlogis", "rt",
  "pexp", "pnorm", "punif", "ppois", "pbinom", "pt", "pchisq", "pf",
  "pweibull", "pgamma", "pnbinom", "pgeom", "phyper", "plogis", "pt",
  "qexp", "qnorm", "qunif", "qpois", "qbinom", "qt", "qchisq", "qf",
  "qweibull", "qgamma", "qnbinom", "qgeom", "qhyper", "qlogis", "qt"
)

## a list of functions that are can be used in the models
attached.functions = new.env()

#' attach a function to be used by the models.
#' @param functions a named list of functions, where the name are the functions used
#' in models
#' @param ... named arguments, each specifies a function
#' @details a NULL value removes a function from the list of attached functions.
#' @export
#' 
attach.function <- function(..., functions=NULL) {
  # the ... must be named functions
  fs = c(list(...), functions)
  ns = names(fs)
  if (is.null(ns) || "" %in% ns)
    stop("functions must be named")
  check = sapply(fs, function(f) is.function(f) || is.null(f))
  if (!all(check))
    stop(paste(ns[!check], collapse=", "), "must be functions")
  for (n in ns)
    attached.functions[[n]] <- fs[[n]]
}

#' A R6 class that is the super class for all numerical simulators
#' 
#' The main method is `simulate`.
#' @docType class
#' @export
Simulator = R6Class(
  "Simulator",
  private = list(
    # a representation of the model
    .model = NULL,
    # model compartments
    compartments = NULL,
    # model parameters
    parameters = NULL,
    # model substitutions
    alias = NULL,
    # attached functions
    attached.functions = NULL,
    # a function call to calculate alias
    f.alias = NULL,
    
    # a private method that builds the code to simulate the
    # model. It should be implemented by all subclasses.
    # Parameter: model the model to simulate. 
    # Return value: the value that will be stored in the .model private field. 
    # This value may then be used in the private method `run` to actually
    # simulate. 
    build = function(model) { 
      NULL 
    },
    
    # format the ODEs as R commands
    format.equation = function(eq) {
      var = eq[[2]]
      if (is.call(var)) {
        if (var[[1]] != "'")
          stop("Invalid equation ", eq)
        var = as.name(paste0(".d.", var[[2]]))
      }
      call("<-", var, eq[[3]])
    },
    
    # create a list of R assignment calls to convert the components of
    # a vector (given by name) into variables with the same name.
    format.var = function(S, name) {
      lapply(S, function(var) {
        call("<-", as.name(var), call("[[", as.name(name), which(var == S)))
      })
    },

    # create a list of R assignment calls to calculate the substitutions
    format.substitution = function() {
      c(
        private$format.var(private$compartments, "y"),
        private$format.var(private$parameters, "parms"),
        lapply(private$alias, private$format.equation)
      )
    },
    
    # This private method actually performs the simulation
    run = function(t, y0, parms, ...) {
      NULL
    },
    
    # run the simulation, and if requested (alias=TRUE), then calculate
    # the substitutions that depend on compartments and append them
    # to the data frame returned from `run`.
    .simulate = function(t, y0, parms, alias, ...) {
      data = private$run(t, y0[private$compartments], parms[private$parameters], ...)
      d = if (!alias) data else
        as.data.frame(t(apply(data, 1, private$f.alias, parms=parms)))
    }
  ),
  
  public = list(
    #' @description constructor
    #' @param model the model to simulate
    #' @details the Simulation class is an abstract class. It's constructor
    #' should not be explicitly called.
    #' 
    #' Functions used by the model must be explicited attached using the attach.function 
    #' function provided by this package.
    initialize = function(model) {
      # model functions must be defined
      fs = model$functions
      if (!is.null(fs)) {
        check = sapply(fs, function(f) f %in% default.functions || !is.null(attached.functions[[f]]))
        if (!all(check)) {
          missing = fs[!check]
          stop("missing functions: ", paste(missing, collapse=", "))
        }
      }
      private$compartments = model$compartments
      private$parameters = model$parameters
      # note that the order of subst guarantees that the evaluation can be carried out
      # sequentially
      subst = model$substitutions
      private$alias = list()
      for (n in names(subst)) {
        s = subst[[n]]
        private$alias[[n]] = list(as.name("<-"), as.name(n), s)
        if (!is.null(attr(s, "compartment")))
          attr(private$alias[[n]], "compartment") = TRUE
      }
      private$.model = private$build(model)
      environment(private$.model) = attached.functions
      col = list(quote(c), quote(row))
      for (n in names(private$alias)) {
        if (!is.null(attr(private$alias[[n]], "compartment")))
          col[[n]] = as.name(n)
      }
      body = call(
        "with", 
        quote(c(as.list(parms), as.list(row))),
        as.call(c(as.name("{"), lapply(private$alias, as.call), as.call(col)))
      )
      private$f.alias = as.function(c(alist(row=, parms=), body))
    },
    
    #' @description simulate the model
    #' @param t a numeric vector of times.
    #' @param y0 a numeric vector of initial conditions.
    #' @param parms a numeric vector of parameter values.
    #' @param vars a character vector specifying the variable names to be returned. 
    #' A variable can be either a compartment or a substitution
    #' @param ... extra arguments to be passed to the `ode` method
    #' substitutions that depend on the compartments are calculated 
    #' and returned
    #' @return a data frame which rows correspond to each time in `t`, 
    #' and columns correspond to the variables specified in vars.
    simulate = function(t, y0, parms=NULL, vars=names(y0), ...) {
      ny = names(y0)
      if (is.null(ny) || any (ny == ""))
        stop("y0 must be named")
      np = names(parms)
      if (is.null(np) || any (np == ""))
        stop("parms must be named")
      y0 = y0[private$compartments]
      missing = which(is.na(y0))
      if (length(missing) > 0)
        stop("missing initial values for ", private$compartments[missing])
      extra = setdiff(private$compartments, ny)
      if (length(extra) > 0)
        stop("extra initial values for ", extra)
      parms = parms[private$parameters]
      missing = which(if (is.list(parms)) sapply(parms, is.null) else is.na(parms))
      if (length(missing) > 0)
        stop("missing parameter values for ", private$parameters[missing])
      extra = setdiff(private$parameters, np)
      if (length(extra) > 0)
        stop("extra parameter values for ", extra)
      alias = any(vars %in% names(private$alias))
      data = private$.simulate(t, y0, parms, alias, ...)
      if (colnames(data)[[1]] == "time")
        vars = c("time", vars)
      data[, vars]
    }
  ),
  
  active = list(
    #' @field model a read-only field returning the  simulation code 
    #' for the model.
    model = function() {
      private$.model
    }
  )
)
