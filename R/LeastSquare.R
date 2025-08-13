#' The calibrator using the least square method.
#' @name LeastSquare
#' @docType class
#' @export
LeastSquare <- R6::R6Class(
  "LeastSquare",
  inherit = Calibrator,
  
  private = list(
    # whether to fit in log scale
    .log = FALSE,

    simulator = function(model) {
      ODE$new(model)
    },
    
    objective = function(pars, formula, fixed, ...) {
      x = private$simulate(pars, formula, fixed, ...)
      if (private$.log) x = log(x)
      if (is.data.frame(x)) {
        sum(sapply(1:ncol(x), function(i) {
          sum((x[,i]-private$.data[, i])^2)
        }))
      } else sum((x-private$.data)^2)
    },
    
    .calibrate = function(guess, formula, fixed, control=NULL, ...) {
      args = list(...)
      if (!is.null(control) && args$method == "Nelder-Mead") {
        control$trace = 1
        control$parscale=guess/10
      }
      optim(guess, private$objective, formula=formula, fixed=fixed, control=control, ...)
    },
    
    interpret = function(result) {
      if (result$convergence == 0) {
        p = result$par
        if (!is.null(result$hessian)) {
          V = solve(result$hessian)
          x = sapply(names(p), function(n) {
            m = p[[n]]
            s = sqrt(V[n, n])
            c(mean=m, "2.5%"=m-2*s, "97.5%"=m+2*s)
          })
          as.data.frame(t(x))
        } else {
          data.frame(mean=p)
        }
      } else {
        print(result)
        stop("Error (", result$convergence, "): ", result$message)
      }
    }
  ),
  
  public = list (
    #' @description initializer
    #' @param model the model to calibrate
    #' @param time either a numeric vector containing the times (including the 
    #' initial time) of the ODE solution that corresponds to the data, or a 
    #' character value giving the name of the column in data that corresponds 
    #' to time.
    #' @param data a data.frame object containign the data for the calibration
    #' @param ... each argument is a formula defining the maps between 
    #' the data columns and the model variables. Please see the details section.
    #' @param log boolean indicating whether to fit in log scale
    #' @param cumulative whether the data is cumulative
    #' @param mapping a list specifying the mapping from data columns to model variables.
    #' @details 
    #' A mapping is a named argument, where name is the
    #' data colummn name, and value corresponds to the model variables (or an 
    #' expression to calculate from the model variables.)
    initialize = function(model, time, data, ..., log=FALSE, cumulative=FALSE, mapping=character()) {
      private$.log = log
      super$initialize(model, time, data, ..., cumulative = cumulative, mapping = mapping)
      if (log) private$.data = base::log(private$.data)
    }
  )
)
