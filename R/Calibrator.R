#' The base R6 class for calibrators
#' 
#' @details A calibrator calibrates the model to a dataset. Like a simulator, 
#' to calibrate the model to a dataset, a calibrator
#' expects the known and unknown initial conditions and the parameter values, 
#' and the data to calibrate to. In addition, it also expects a mapping from 
#' model solutions to observation variables.
#' @name Calibrator
#' @docType class
#' @export

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
    # the return value from the underlying algorithm
    .details = NULL,
    
    ## create a simulator needed for calibration. Must be implemented by subclasses.
    simulator = function(model) {
      NULL
    },
    
    ## the function that simulates the model given the fit info.
    simulate = function(pars, formula, fixed, ...) {
      # evaluate values given by the formula
      all = c(as.list(pars), fixed)
      for (n in names(formula)) {
        f = formula[[n]]
        all[[n]] = eval(f$expr, envir = all)
      }
      ic = unlist(all[private$.model$compartments])
      par = all[private$.model$parameters]
      data = private$.simulator$simulate(private$.time, ic, par, vars=unlist(private$.mapping), ...)
      if (private$.cumulative) {
        if (ncol(data) == 2) diff(data[,2]) else
          as.data.frame(lapply(data[, -1], diff))
      } else data[,-1]
    },
    
    ## interpret the results of calibration
    interpret = function(results) {
      NULL
    },
    
    ## The actual calibration is done in this function.
    ## This method Must be implemented by subclasses.
    ## @param guess the initial guess for the parameters to be fitted
    ## @param formula gives the parameter values that should be calculated
    ## using this formula
    ## @param fixed the parameter values that are given
    ## @return the fitting results, which will be passed to 
    ## the interpret method.
    .calibrate = function(guess, formula, fixed, ...) {
      NULL
    },
    
    # split the parameters or initial conditions by type
    # here x is the list or vector of parameters or initial conditions,
    # mode specifies what values to split. If x does not 
    # contain a named value, it will be fitted
    split = function(x, mode=c("initial.value", "parameter")) {
      mode = match.arg(mode, c("initial.value", "parameter"))
      if (mode == "initial.value") {
        set = private$.model$compartments
      } else if (mode == "parameter") {
        set = private$.model$parameters
      } else stop("invalid mode: ", mode)
      ns = names(x)
      if (is.null(ns) || any(ns == ""))
        stop(mode, "must be named")
      extra = setdiff(ns, set)
      if (length(extra) > 0)
        stop(paste(extra, collapse = ", "), " not defined in the model")
      idx.values = sapply(x, is.numeric)
      idx.formula = sapply(x, function(y) is(y, "Expression"))
      x.values = as.list(x[idx.values])
      x.formula = x[idx.formula]
      x.fit = setdiff(set, ns[idx.values | idx.formula])
      list(value = x.values, formula = x.formula, fit = x.fit)
    },
    
    # the calibrate method uses this function to infer the parameters that 
    # needs to be fitted, those that need to be calculated from a formula, 
    # and those that are fixed.
    # this function returns a list with names formula, fit and fixed
    # to be fitted, in addition to the remaining parameters in ...
    fit.info = function(initial.values, parms, guess, ...) {
      ic = private$split(initial.values, mode="initial.value")
      p = private$split(parms, mode="parameter")
      formula = c(ic$formula, p$formula)
      if (length(formula) > 0) formula = formula[order(formula)]
      # find all parameters from the formula
      parms = Reduce(
        function(extra, f) union(f$parms, extra), 
        formula, 
        init=character())
      extra = setdiff(parms, c(private$.model$compartments, private$.model$parameters))
      fit = c(extra, ic$fit, p$fit)
      if (length(fit) == 0)
        stop("no initial values or parameters to fit")
      if (!is.numeric(guess) || any(is.na(guess)))
        stop("invalid initial values")
      ng = names(guess)
      if (is.null(ng) || any(ng == ""))
        stop("guess must be a named list")
      extra = setdiff(ng, fit)
      if (length(extra) == 1) {
        stop("guess contains extra parameter: ", extra)
      } else if (length(extra) > 1) {
        stop("guess contains extra parameters: ", paste(extra, collapse=", "))
      }
      missed = setdiff(fit, ng)
      if (length(missed) == 1) {
        stop("guess is missing parameter: ", missed)
      } else if (length(missed) > 1) {
        stop("guess is missing parameters: ", paste(missed, collapse=", "))
      }
      fixed = c(ic$value, p$value)
      c(
        list(
          guess = guess,
          formula = formula,
          fixed = fixed
        ),
        list(...)
      )
    }
  ),
  
  public = list(
    #' @description initializer
    #' @param model the model to calibrate
    #' @param time either a numeric vector containing the times (including the 
    #' initial time) of the ODE solution that corresponds to the data, or a 
    #' character value giving the name of the column in data that corresponds 
    #' to time.
    #' @param data a data.frame object containing the data for the calibration
    #' @param ... each argument is a formula defining the maps between 
    #' the data columns and the model variables. Please see the details section.
    #' @param cumulative whether the data is cumulative
    #' @param mapping a named list specifying the mapping from data columns
    #' to quantities in the model.
    #' @details 
    #' A mapping is a named argument, where name is the
    #' data colummn name, and value corresponds to the model variables (or an 
    #' expression to calculate from the model variables.)
    initialize = function(model, time, data, ..., cumulative=FALSE, mapping=character()) {
      m = model$clone(deep=TRUE)
      if (!is.data.frame(data))
        stop("data must be a data.frame object")
      extra = setdiff(names(mapping), colnames(data))
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
      private$.data = data[, names(private$.mapping)]
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
      info = private$fit.info(initial.values, parms, guess, ...)
      private$.details = do.call(private$.calibrate, info)
      private$interpret(private$.details)
    }
  ),
  
  active = list(
    #' @field parameters the names of all parameters
    parameters = function() {
      private$.model$parameters
    },
    
    #' @field details the details of the fitting output
    details = function() {
      private$.details
    }
  )
)
