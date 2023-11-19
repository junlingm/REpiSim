MLE <- R6Class(
  "MLE",
  inherit = Optimizer,
  
  private = list(
    .likelihood = NULL,
    
    parameters = function() {
      c(private$.model$parameters, private$.likelihood$par)
    },
    
    objective = function(pars, initial.values, parms) {
      pars.l = pars[private$.likelihood$par]
      pars.l = pars.l[!is.na(pars.l)]
      logL = private$.likelihood$logL
      x = private$simulate(pars, initial.values, parms)
      if (is.data.frame(x)) {
        -sum(sapply(1:ncol(x), function(i) {
          logL(private$.data[,i], x[,i], pars.l)
        }))
      } else -logL(private$.data[,1], x, pars.l)
    },
    
    optimizer = function(pars, initial.values, parms, neglogL, ...) {
      args = names(pars)
      arglist = list(as.name("c"))
      for (arg in args) arglist[[arg]] = as.name(arg)
      body = call("{",
        call("<-", as.name("pars"), as.call(arglist)),
        quote(private$objective(pars, initial.values, parms))
      )
      neglogL = as.function(c(
        as.list(pars),
        body
      ))
      e = new.env()
      e$private = private
      environment(neglogL)=e
      mle2(neglogL, as.list(pars), ...)
    },
    
    interpret = function(result) {
      details = result@details
      if (details$convergence == 0) {
        ci = confint(result)
        if ("mle2" %in% class(ci)) {
          x = coef(ci)
          n = names(coef(result))
          error = paste(sapply(1:length(x), function(i) {
            paste0(n[[i]], "=", x[[i]])
          }), collapse=", ")
          stop("a better result is found: ", error)
        }
        as.data.frame(cbind(mean=coef(result), ci))
      } else stop("Error (", details$convergence, "): ", details$message)
    },
    
    .calibrate = function(pars, ...) {
      miss = setdiff(private$.likelihood$par, pars)
      if (length(miss) > 0)
        stop("missing distrbution parameter", if (length(miss)==1) "" else "s", 
             ":", miss)
      super$.calibrate(pars, ...)
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
    #' @details 
    #' A mapping is a named argument, where name is the
    #' data colummn name, and value corresponds to the model variables (or an 
    #' expression to calculate from the model variables.)
    initialize = function(model, time, data, likelihood, ..., cumulative=FALSE, mapping=character()) {
      library(bbmle)
      if (!"Distribution" %in% class(likelihood)) 
        stop("likelihood must be a Distribution object")
      private$.likelihood = likelihood$likelihood()
      super$initialize(model, time, data, ..., cumulative = cumulative, mapping = mapping)
    }
  )
)
