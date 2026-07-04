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
    .last = NULL, # th e last fit result
    
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
      caller = list(...)
      if (is.null(caller$skip.hessian))
        caller$skip.hessian=!private$.CI
      res = do.call(bbmle::mle2, c(list(neglogL, start=as.list(guess)), caller))
      if (res@details$convergence == 10 && (is.null(private$.last) || any(coef(res)!=private$.last))) {
        private$.last <- coef(res)
        cat("The likelihood profiling did not converge. continue...\n")
        print(coef(res))
        private$.calibrate(coef(res), formula, fixed, ...)
      } else res
    },
    
    interpret = function(result) {
      details = result@details
      if (details$convergence == 0) {
        estimate = bbmle::coef(result)
        if (private$.CI) {
          ci = bbmle::confint(result)
          if ("mle2" %in% class(ci)) {
            x = bbmle::coef(ci)
            n = names(estimate)
            error = paste(sapply(1:length(x), function(i) {
              paste0(n[[i]], "=", x[[i]])
            }), collapse=", ")
            stop("a better result is found: ", error)
          }
        } else ci = NULL
        as.data.frame(cbind(mean=estimate, ci))
      } else {
        .last.fit <<- result
        cat("stopped at\n")
        print(bbmle::coef(result))
        stop("Error (", details$convergence, "): ", details$message)
      }
    },
    
    fit.info = function(initial.values, parms, guess, ...) {
      # split private$.likelihood$par from parms
      parms = as.list(parms)
      v = parms[private$.par.likelihood]
      for (i in private$.par.likelihood) {
        parms[[i]] = NULL
      }
      
      likelihood.fit = character(0)
      likelihood.fixed = list()
      likelihood.formula = list()
      
      for (i in private$.par.likelihood) {
        if (calibrator.is.missing(v[[i]])) {
          likelihood.fit = c(likelihood.fit, i)
        } else if (calibrator.is.formula(v[[i]])) {
          likelihood.formula[[i]] = calibrator.as.expression(v[[i]])
        } else if (is.numeric(v[[i]])) {
          likelihood.fixed[[i]] = v[[i]]
        } else {
          stop("Invalid value for likelihood parameter ", i, ": ", v[[i]])
        }
      }
      
      likelihood.extra = if (length(likelihood.formula) > 0) {
        needed = Reduce(function(extra, f) union(f$parms, extra), likelihood.formula, init = character())
        setdiff(needed, c(private$.model$compartments, private$.model$parameters, private$.par.likelihood))
      } else character(0)
      likelihood.fit = c(likelihood.extra, likelihood.fit)
      
      model.guess = guess[setdiff(names(guess), c(private$.par.likelihood, likelihood.extra))]
      info = super$fit.info(initial.values, parms, model.guess, ...)
      
      unexpected = intersect(names(guess), setdiff(private$.par.likelihood, likelihood.fit))
      if (length(unexpected) == 1) stop("guess contains extra parameter: ", unexpected)
      if (length(unexpected) > 1) stop("guess contains extra parameters: ", paste(unexpected, collapse = ", "))
      
      missed = setdiff(likelihood.fit, names(guess))
      if (length(missed) == 1) stop("guess is missing parameter: ", missed)
      if (length(missed) > 1) stop("guess is missing parameters: ", paste(missed, collapse = ", "))
      
      info$guess = guess
      info$formula = c(info$formula, likelihood.formula)
      info$fixed = c(info$fixed, likelihood.fixed)
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
      if (!requireNamespace("bbmle", quietly = TRUE))
        stop("bbmle package is required for MLE calibration")
      if (!"Distribution" %in% class(likelihood) || is.null(likelihood$log.likelihood))
        stop("likelihood must be a Distribution object")
      private$.likelihood = likelihood
      private$.par.likelihood = likelihood$parameters
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
