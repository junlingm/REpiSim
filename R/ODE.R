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
    },

    #' @description calculate the equilibrium of the model
    #' 
    #' Calculate the equilibrium of the the solution with a given initial
    #' condition.
    #' 
    #' @param time a numeric vector giving times.
    #' @param y0 a named numeric vector giving initial condition.
    #' @param parms a named numeric vector giving the parameter values.
    #' @param max.iter the maximum iterations to search for the equilibrium
    #' @param error a numeric value giving the 2-norm of the tolerance 
    #' @param vars a character vector giving the compartments and substitutions
    #' which equilibrium value should be returned.
    #' @param ... extra parameters to be passed to the `deSolve::ode` function.
    #' @return a named numeric vector giving the equilibrium value, or NULL
    #' if convergence is not achieved.
    #' @details the method repetitively call the simulate method to solve
    #' the system using the given `time`, `y0`, and `parms`. The 2-norm of the
    #' difference between y0 and the final state is calculated. If it is less then 
    #' `error` then the final state is returned as the equilibrium. Otherwise,
    #' the initial state is reset to be the final state, and the process is
    #' repeated for at most `max.iter` times.
    #' @examples 
    #'     #' # an SIR model
    #' SIR = Model$new(title="SIR")
    #' SIR$compartment(S ~ -beta*S*I/N)$
    #'   compartment(I ~ beta*S*I/N - gamma*I)$
    #'   compartment(R ~ gamma*I)$
    #'   where(N = S + I + R)
    #' ode = ODE$new(SIR)
    #' ode$equilibrium(
    #'   time=0:40, 
    #'   y0=c(S=10000, I=1, R=0), 
    #'   parms=c(beta=0.4, gamma=0.2), 
    #'   alias = FALSE) # N is not returned.
    equilibrium = function(time, y0, parms, max.iter = 100, error = 1e-6, vars = names(y0), ...) {
      n = length(time)
      m = length(range)
      y = NULL
      alias = any(sapply(vars, function(x) !is.null(private$alias[[x]])))
      for (i in 1:max.iter) {
        sol = private$.simulate(time, y0, parms, alias = alias, ...)
        k = ncol(sol)-1
        remove = if (k > length(y0)) (length(y0)+1):k else c()
        y = unlist(sol[n, -1])
        y1 = if (length(remove) > 0) y[-remove] else y
        if (sqrt(sum((y0[vars]-y1[vars])^2, na.rm = TRUE)) < error)
          break
        y0 = y1
      }
      y[vars]
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
