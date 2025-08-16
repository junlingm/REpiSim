#' A maximum likelihood calibrator using the bbmle package
#' @name MLE
#' @docType class
#' @export
MLE <- R6Class(
  "MLE",
  inherit = Calibrator,
  
  private = list(
    .likelihood = NULL,
    .par.likelihood = NULL,
    .CI = TRUE,
    
    objective = function(pars, formula, fixed, ...) {
      all = c(as.list(pars), fixed)
      for (n in names(formula)) {
        f = formula[[n]]
        all[[n]] = eval(f$expr, envir = all)
      }
      pars.l = all[private$.par.likelihood]
      x = private$simulate(all, NULL, NULL, ...)
      if (is.data.frame(x)) {
        -sum(sapply(1:ncol(x), function(i) {
          do.call(
            private$.likelihood$log.likelihood, 
            c(list(private$.data[,i], x[,i]), pars.l)
          )
        }))
      } else -do.call(
        private$.likelihood$log.likelihood,
        c(list(private$.data, x), pars.l)
      )
    },
    
    simulator = function(model) {
      ODE$new(model)
    },
    
    .calibrate = function(guess, formula, fixed, ...) {
      args = names(guess)
      arglist = list(as.name("c"))
      for (arg in args) arglist[[arg]] = as.name(arg)
      body = call("{",
        call("<-", as.name("pars"), as.call(arglist)),
        quote(private$objective(pars, formula, fixed))
      )
      neglogL = as.function(c(
        as.list(guess),
        body
      ))
      e = new.env()
      e$fixed = fixed
      e$formula = formula
      e$private = private
      environment(neglogL)=e
      mle2(neglogL, as.list(guess), ...)
    },
    
    interpret = function(result) {
      details = result@details
      if (details$convergence == 0) {
        if (private$.CI) {
          ci = confint(result)
          if ("mle2" %in% class(ci)) {
            x = coef(ci)
            n = names(coef(result))
            error = paste(sapply(1:length(x), function(i) {
              paste0(n[[i]], "=", x[[i]])
            }), collapse=", ")
            stop("a better result is found: ", error)
          }
        } else ci = NULL
        as.data.frame(cbind(mean=coef(result), ci))
      } else stop("Error (", details$convergence, "): ", details$message)
    },
    
    fit.info = function(initial.values, parms, ...) {
      # split private$.likelihood$par from parms
      parms = as.list(parms)
      v = parms[private$.par.likelihood]
      for (i in private$.par.likelihood) {
        parms[[i]] = NULL
      }
      info = super$fit.info(initial.values, parms, ...)
      for (i in private$.par.likelihood) {
        if (is.na(v[[i]]) || is.null(v[[i]])) {
          info$fit = c(info$fit, i)
        } else if (is(v[[i]], "Expression")) {
          info$formula = c(info$formula, v[i])
        } else if (is.numeric(v[[i]])) {
          info$fixed[[i]] = v[[i]]
        } else {
          stop("Invalid value for likelihood parameter ", i, ": ", v[[i]])
        }
      }
      info
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
    #' @param likelihood a Lieklihood object specify the type of likelihood
    #' @param ... each argument is a formula defining the maps between 
    #' the data columns and the model variables. Please see the details section.
    #' @param cumulative whether the data is cumulative
    #' @param mapping a named list specifying the mapping from data columns
    #' @param CI a boolean indicating whether to calculate the CI using likelihood profiling
    #' @details 
    #' A mapping is a named argument, where name is the
    #' data colummn name, and value corresponds to the model variables (or an 
    #' expression to calculate from the model variables.)
    initialize = function(model, time, data, likelihood, ..., cumulative=FALSE, mapping=character(), CI=TRUE) {
      if (!require(bbmle))
        stop("bbmle package is required for MLE calibration")
      if (!"Distribution" %in% class(likelihood) || is.null(likelihood$log.likelihood))
        stop("likelihood must be a Distribution object")
      private$.likelihood = likelihood
      l = formals(likelihood$log.likelihood)
      if (length(l) >= 3) private$.par.likelihood = names(l[3:length(l)])
      private$.CI = CI
      super$initialize(model, time, data, ..., cumulative = cumulative, mapping = mapping)
    }
  ),
  
  active = list(
    #' @field parameters the names of all parameters
    parameters = function() {
      c(private$.model$parameters, private$.likelihood$par)
    },
    
    #' @field CI whether to calculate the CI using likelihood profiling
    CI = function() {
      private$.CI
    }
  )
)
