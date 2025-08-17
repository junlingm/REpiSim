#' An MCMC calibrator that uses the mcmc::metrop function
#' @name Metrop
#' @docType class
#' @export
Metrop <- R6::R6Class(
  "Metrop",
  inherit = Baysian,
  
  private = list(
    .likelihood = NULL,
    
    objective = function(pars, formula, fixed, priors, ...) {
      names(pars) = names(priors)
      all = c(as.list(pars), fixed)
      lp = sum(sapply(names(priors), function(n) priors[[n]]$log.density(all[[n]]) ))
      if (is.infinite(lp) || is.nan(lp)) {
        return(-Inf)
      }
      for (n in names(formula)) {
        f = formula[[n]]
        all[[n]] = eval(f$expr, envir = all)
      }
      pars.l = all[private$.par.likelihood]
      x = private$simulate(all, NULL, NULL, ...)
      ll = if (is.data.frame(x)) {
        sum(sapply(1:ncol(x), function(i) {
          do.call(
            private$.likelihood$log.likelihood, 
            c(list(private$.data[,i], x[,i], pars.l))
          )
        }))
      } else do.call(
        private$.likelihood$log.likelihood,
        c(list(private$.data, x), pars.l)
      )
      ll + lp
    },
    
    .calibrate = function(guess, formula, fixed, priors, ...) {
      x = metrop(private$objective, guess, ..., 
                 formula=formula,
                 fixed=fixed,
                 priors = priors)
      colnames(x$batch) = names(guess)
      x
    },
    
    interpret = function(result) {
      s = summary(as.mcmc(result$batch), quantiles=c(0.025, 0.975))
      as.data.frame(cbind(
        mean = s$statistics[,"Mean"],
        s$quantiles
      ))
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
    #' @param likelihood a Distribution object specify the type of likelihood
    #' @param ... each argument is a formula defining the maps between 
    #' the data columns and the model variables. Please see the details section.
    #' @param cumulative whether the data is cumulative
    #' @param mapping a list specifying the mapping from data columns to model variables.
    #' @details 
    #' A mapping is a named argument, where name is the
    #' data colummn name, and value corresponds to the model variables (or an 
    #' expression to calculate from the model variables.)
    initialize = function(model, time, data, likelihood, ..., cumulative=FALSE, mapping=character()) {
      library(mcmc)
      library(coda)
      if (!"Distribution" %in% class(likelihood)) 
        stop("likelihood must be a Distribution object")
      private$.likelihood = likelihood
      super$initialize(model, time, data, ..., cumulative = cumulative, mapping = mapping)
    },
    
    #' continue from the previous run
    #' @param ... extra arguments to be passed to the metrop function
    continue = function(...) {
      if (is.null(private$.details))
        stop("was not run before")
      guess = private$.details$batch[nrow(private$.details$batch), ]
      private$.details = private$metrop(
        private$.details, 
        guess=guess,
        initial.values=private$.details$initial.values, 
        parms=private$.details$parms, ...)
      private$interpret(private$.details)
    }
  ),
  
  active = list(
    #' @field parameters the names of all parameters
    parameters = function() {
      c(private$.model$parameters, private$.likelihood$par)
    },
    
    #' @field samples the samples of the posterior distribution
    samples = function() {
      as.mcmc(private$.details$batch)
    }
  )
)
