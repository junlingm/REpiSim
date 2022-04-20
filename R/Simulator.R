#' A R6 class that is the super class for all numerical simulators
#' 
#' The main method is `simulate`.
#' @docType class
#' @exportClass Simulator
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
    
    # This private method actually performs the simulation
    run = function(t, y0, parms, ...) {
      NULL
    },
    
    # run the simulation, and if requested (alias=TRUE), then calculate
    # the substitutions that depend on compartments and append them
    # to the data frame returned from `run`.
    .simulate = function(t, y0, parms, alias, ...) {
      data = private$run(t, y0[private$compartments], parms[private$parameters], ...)
      if (!alias) data else
        as.data.frame(t(apply(data, 1, private$f.alias, parms=parms)))
    }
  ),
  
  public = list(
    #' @description constructor
    #' @param model the model to simulate
    #' @details the Simulation class is an abstract class. It's constructor
    #' should not be explicitly called.
    initialize = function(model) {
      private$compartments = model$compartments
      private$parameters = model$parameters
      private$.model = private$build(model)
      subst = model$substitutions
      private$alias = list()
      for (n in names(subst)) {
        s = subst[[n]]
        private$alias[[n]] = private$format.equation(list("=", as.name(n), s))
        if (!is.null(attr(s, "compartment")))
          attr(private$alias[[n]], "compartment") = TRUE
      }
      col = list(quote(c), quote(row))
      for (n in names(private$alias)) {
        if (!is.null(attr(private$alias[[n]], "compartment")))
          col[[n]] = as.name(n)
      }
      body = call(
        "with", 
        quote(c(as.list(parms), as.list(row))),
        as.call(c(as.name("{"), private$alias, as.call(col)))
      )
      private$f.alias = as.function(c(alist(row=, parms=), body))
    },
    
    #' @description simulate the model
    #' @param t a numeric vector of times.
    #' @param y0 a numeric vector of initial conditions.
    #' @param parms a numeric vector of parameter values.
    #' @param alias a logical value. If TRUE, the values of the 
    #' @param ... extra arguments to be passed to the `ode` method
    #' substitutions that depend on the compartments are calculated 
    #' and returned
    #' @return a data frame which rows correspond to each time in `t`, 
    #' and columns correspond to the value of each compartment (or
    #' substitution if alias=TRUE).
    simulate = function(t, y0, parms=NULL, alias=TRUE, ...) {
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
      private$.simulate(t, y0, parms, alias, ...)
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
