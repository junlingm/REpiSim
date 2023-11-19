#' An R6 class to calibrate a model to data
#' 
#' @details Like a simulator, to calibrate the model to a dataset, a calibrator
#' expects the known and unknown initial conditions and the parameter values, 
#' and the data to calibrate to. In addition, it also expects a mapping from 
#' model solutions to observation variables.
Calibrator <- R6::R6Class(
  "Calibrator",
  private = list(
    # the simulator to simulator the model udrign calibration
    .simulator = NULL,
    # the data frame to be fitted to
    .data = NULL,
    # the model fitted to data
    .model = NULL,
    # whether the solutions of the model to be matched to data is cumulative
    .cumulative = FALSE,
    # the mapping from model variables to data
    .mapping = NULL,
    # the time for simulation
    .time = NULL,
    
    ## create a simulator needed for calibration. Must be implemented by subclasses.
    simulator = function(model) {
      NULL
    },
    
    ## check parameters. return the names of parameters to be fitted
    parameters = function() {
      private$.model$parameters
    },

    simulate = function(fit, ic, parms) {
      ic.filled = ic$value
      ic.filled[ic$fit] = fit[ic$fit]
      pars.filled = parms$value
      pars.filled[parms$fit] = fit[parms$fit]
      data = private$.simulator$simulate(private$.time, ic.filled, pars.filled, vars=private$.mapping)
      if (private$.cumulative) {
        l = list()
        n = names(data)
        for (col in 2:nrow(data)) {
          l[[n[[col]]]] = diff(data[,col])
        }
        as.data.frame(l)
      } else data[,-1]
    },
    
    ## interpret the results of calibration
    interpret = function(results) {
      NULL
    },
    
    ## The actual calibration is done in this function.
    ## This method Must be implemented by subclasses.
    .calibrate = function(pars, intial.values, parms, ...) {
      NULL
    }
  ),
  
  public = list(
    #' @description initializer
    #' @param model the model to calibrate
    #' @param time either a numeric vector containing the times (including the 
    #' initial time) of the ODE solution that corresponds to the data, or a 
    #' character value giving the name of the column in data that corresponds 
    #' to time.
    #' @param data a data.frame object containign the data for the calibration
    #' @param ... each argument is a formula defining the maps between 
    #' the data columns and the model variables. Please see the details section.
    #' @details 
    #' A mapping is a named argument, where name is the
    #' data colummn name, and value corresponds to the model variables (or an 
    #' expression to calculate from the model variables.)
    initialize = function(model, time, data, ..., cumulative=FALSE, mapping=character()) {
      m = model$clone(deep=TRUE)
      if (!is.data.frame(data))
        stop("data must be a data.frame object")
      extra = setdiff(names(mapping), names(data))
      if (length(extra) > 0)
        stop("data column", if(length(extra)==1) "" else "s", 
             " does not exist: ", paste(extra, collapse=", "))
      private$.mapping = mapping
      private$.cumulative = cumulative
      if (is.numeric(time)) {
        if (cumulative) {
          if (length(time) - 1 != nrow(data))
            stop("the time should have one more point than the length of the data, to calculate the difference from teh cumulative curves")
        } else if (length(time) != nrow(data))
          stop("time should have the same length as data")
        private$.time = time
        l = list()
        private$.data= do.call(data.frame, l)
      } else if (!is.character(time) || is.null(data[[time]]))
        stop("time should be a column name in data indicating which column is time")
      else if (cumulative)
        stop("the time should have one more point than the length of the data, to calculate the difference from teh cumulative curves")
      else private$.time = data[[time]]
      args = as.list(substitute(list(...))[-1])
      if (length(args) == 0 && length(mapping) == 0)
        stop("no mapping from data to model variables")
      for (a in args) {
        if (as.character(a[[1]]) != "~")
          stop("invalid mapping: ", deparse(a))
        if ((!is.name(a[[2]]) && !is.character(a[[2]])) || is.null(data[[a[[2]]]])) 
          stop("not a data column: ", deparse(a[[2]]))
        col = as.character(a[[2]])
        if (is.name(a[[3]])) {
          var = as.character(a[[3]])
          if (!var %in% m$compartments && is.null(m$substitutions[[var]]))
            stop(var, "is not a model variable")
        } else {
          if (!is.call(a[[3]]))
            stop(deparse(a[[3]]), " is not a valid specification")
          var = if (!col %in% m$compartments && is.null(m$substitutions[[col]])) {
            col
          } else paste0(".observation.", col)
          l = list()
          l[[var]] = a[[3]]
          m$where(pairs = l)
        }
        if (!is.na(private$.mapping[col]))
          stop("redefining mapping for data column ", col)
        private$.mapping[col] = var
      }
      private$.model = m
      private$.simulator = private$simulator(m)
      private$.data = data[names(private$.mapping)]
    },
    
    #' Calibrate the model to data
    #' 
    #' @param initial.values the initial values for the model. The parameters 
    #' that need to be estimate should be NA, those that do not need to be
    #' estimated must contain a finite value.
    #' @param parms the parameter values of the model. The parameters 
    #' that need to be estimate should be NA, those that do not need to be
    #' estimated must contain a finite value.
    #' @param ... extra arguments to be passed to calibrators
    calibrate = function(initial.values, parms, ...) {
      initial.values = initial.values[!is.na(initial.values)]
      n = names(initial.values)
      if (is.null(n) || any(n=="")) 
        stop("initial values must be named")
      extra = setdiff(n, private$.model$compartments)
      if (length(extra) > 0)
        stop("variable", if(length(extra)==1) "" else "s", 
             " not defined in model: ", paste(extra, collapse=", "))
      pars.ic = setdiff(private$.model$compartments, n)
      parms = parms[!is.na(parms)]
      n = names(parms)
      if (is.null(n) || any(n=="")) 
        stop("parameter values must be named")
      extra = setdiff(n, private$parameters())
      if (length(extra) > 0)
        stop("parameter", if(length(extra)==1) "" else "s", 
             " not defined in model: ", paste(extra, collapse=", "))
      pars.parms = setdiff(private$parameters(), n)
      model.parms = intersect(n, private$.model$parameters)
      pars = c(pars.ic, pars.parms)
      x = private$.calibrate(pars, 
                             list(value=initial.values, fit=pars.ic),
                             list(value=parms[model.parms], 
                                  fit=intersect(pars.parms, private$.model$parameters)), 
                             ...)
      private$interpret(x)
    }
  )
)
