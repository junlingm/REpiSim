#' A R6 class for numerically solving the ODE models
#' 
#' ODE is a subclass of `Simulator`. This class is only available if the `deSolve` 
#' package is installed.
#' @name ODE
#' @docType class
#' @export
ODE = R6Class(
  "ODE",
  inherit = Simulator,
  private = list(
    build = function(model) {
      l = alist()
      v = list()
      system = model$equations
      der = lapply(names(system$equations), function(var) as.name(paste0(".d.", var)))
      l = c(
        as.name("{"),
        private$format.substitution(),
        lapply(system$equations, private$format.equation),
        call("list", as.call(c(list(as.name("c")), der)))
      )
      as.function(c(alist(t=, y=, parms=), as.call(l)))
    },
    
    run = function(t, y0, parms, ...) {
      as.data.frame(ode(y0, t, self$model, parms, ...))
    }
  ),
  
  public = list(
    #' @description constructor
    #' @param model the model to simulate. It must be an object of a subclass
    #' of `Model`.
    initialize = function(model) {
      if (!require(deSolve, quietly = TRUE))
        stop("package deSolve is not installed. The ODE class is not available")
      super$initialize(model)
    }
  ),
  
  active = list(
    #' @field model a read-only field giving the R function for the ODE equations
    #' to be passed to deSolve::ode as the `func` argument.
    model = function() {
      private$.model
    }
  )
)

#' An R6 class to calculate the equilibrium of a model
#' 
#' Calculate the equilibrium of the model with a given initial condition
#' 
#' @param model a Model object.
#' @param max.iter the maximum iterations to search for the equilibrium
#' @param vary the name of a parameter to value, a character
#' @param range the range of the varying parameter, a numerical vector
#' @param error a numeric value giving the 2-norm of the tolerance
#' @param vars a character vector giving the compartments and substitutions
#' which equilibrium value should be returned.
#' @examples 
#' ## an SIR model
#' SIR = Model$new(title="SIR")
#' SIR$compartment(S ~ -beta*S*I/N)$
#'   compartment(I ~ beta*S*I/N - gamma*I)$
#'   compartment(R ~ gamma*I)$
#'   where(N = S + I + R)
#' eq = Equilibrium$new(SIR, vary="beta", range=seq(0, 0.4, by=0.01))
#' ode$simulate(
#'   time=0:40, 
#'   y0=c(S=10000, I=1, R=0), 
#'   parms=c(beta=0.4, gamma=0.2), 
#'   alias = FALSE) # N is not returned.
#' @name Equilibrium
#' @docType class
#' @export
Equilibrium = R6Class(
  "Equilibrium",
  inherit = ODE,
  private = list(
    vary = NULL,
    range = NULL,
    max.iter = NULL,
    error = NULL,
    
    eq = function(time, y0, parms, ...) {
      n = length(time)
      y = NULL
      for (i in 1:private$max.iter) {
        sol = super$run(time, y0, parms, ...)
        y = unlist(sol[nrow(sol), -1])
        e = sqrt(sum((y0-y)^2))
        if (!is.na(e) && e < private$error)
          return(sol[nrow(sol), -1])
        y0 = y
      }
      y0[1:length(y0)] = NA
      as.data.frame(as.list(y0))
    },
    
    run = function(time, y0, parms, ...) {
      data = data.frame()
      if (is.null(private$vary)) {
        data = rbind(data, private$eq(time, y0, parms))
      } else {
        for (v in private$range) {
          parms[[private$vary]] = v
          data = rbind(data, private$eq(time, y0, parms))
        }
      }
      data
    }
  ),
  public = list(
    #' @description constructor
    #' @param model the model to simulate. It must be an object of a subclass
    #' of `Model`.
    initialize = function(model, vary=NULL, range=NULL, max.iter=100, error = 1e-6) {
      super$initialize(model)
      private$vary = vary
      private$range = range
      if (!is.null(vary)) {
        if (!is.character(vary) && ! vary %in% model$parameters) 
          stop("vary must be a parameter name")
        if (!is.numeric(range)) 
          stop("range must be a numeric vector")
      }
      private$max.iter = max.iter
      private$error = error
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
      data = super$simulate(t, y0, parms, ...)
      if (!is.null(private$vary)) {
        col = list()
        col[[private$vary]] = private$range
        data = cbind(col, data)
      }
      data
    }
  )
)
    