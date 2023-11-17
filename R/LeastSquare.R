LeastSquare <- R6Class(
  "LeastSquare",
  inherit = Calibrator,
  
  private = list(
    # whether to fit in log scale
    .log = FALSE,
    
    simulator = function(model) {
      ODE$new(model)
    },
    
    lsq = function(pars, initial.values, parms) {
      x = private$simulate(pars, initial.values, parms)
      if (private$.log) x = log(x)
      if (is.data.frame(x)) {
        sum(sapply(1:ncol(x), function(i) {
          sum((x[,i]-private$.data[,i])^2)
        }))
      } else sum((x-private$.data[1])^2)
    },
    
    .calibrate = function(pars, initial.values, parms, guess, ...) {
      if (!is.numeric(guess) || any(is.na(guess)))
        stop("invalid initial values")
      ng = names(guess)
      if (is.null(ng) || any(ng=="")) 
        stop("initial guesses must be named")
      extra = setdiff(pars, ng)
      if (length(extra) > 0)
        stop("variable", if(length(extra)==1) "" else "s", 
             " not defined in model: ", paste(extra, collapse=", "))
      miss = setdiff(ng, pars)
      if (length(miss) > 0)
        stop("missing initial guess", if(length(miss)==1) "" else "es", ": ", 
             paste(miss, collapse=", "))
      optim(guess[pars], private$lsq, hessian=TRUE, 
            initial.values=initial.values, parms=parms)
    },
    
    interpret = function(result) {
      if (result$convergence == 0) {
        V = solve(result$hessian)
        p = result$par
        x = sapply(names(p), function(n) {
          m = p[[n]]
          s = sqrt(V[n, n])
          c(mean=m, "2.5%"=m-2*s, "97.5%"=m+2*s)
        })
        as.data.frame(t(x))
      } else stop("Error (", result$convergence, "): ", result$message)
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
    #' @details 
    #' A mapping is a named argument, where name is the
    #' data colummn name, and value corresponds to the model variables (or an 
    #' expression to calculate from the model variables.)
    initialize = function(model, time, data, ..., log=FALSE, cumulative=FALSE, mapping=character()) {
      private$.log = log
      super$initialize(model, time, data, ..., cumulative = cumulative, mapping = mapping)
      if (log) private$.data = base::log(private$.data)
    },
    
    #' Calibrate the model to data
    #' 
    #' @param initial.values the initial values for the model. The parameters 
    #' that need to be estimate should be NA, those that do not need to be
    #' estimated must contain a finite value.
    #' @param parms the parameter values of the model. The parameters 
    #' that need to be estimate should be NA, those that do not need to be
    #' estimated must contain a finite value.
    #' @param guess the initial guess of the parameters to be fitted
    #' @param ... extra arguments to be passed to calibrators
    calibrate = function(initial.values, parms, guess, ...) {
      super$calibrate(initial.values, parms, guess, ...)
    }
  )
)
