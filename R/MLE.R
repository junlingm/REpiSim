likelihood = function(likelihood, data, simulate, pars, initial.values, parms) {
  pars.l = pars[likelihood$par]
  pars.l = pars.l[!is.na(pars.l)]
  if (length(pars.l) == 0) pars.l = parms$value[likelihood$par]
  logL = likelihood$logL
  x = simulate(pars, initial.values, parms)
  if (is.data.frame(x)) {
    -sum(sapply(1:ncol(x), function(i) {
      logL(data[,i], x[,i], pars.l)
    }))
  } else -logL(data[,1], x, pars.l)
}

#' A maximum likelihood calibrator using the bbmle package
MLE <- R6Class(
  "MLE",
  inherit = Optimizer,
  
  private = list(
    .likelihood = NULL,
    
    objective = function(pars, initial.values, parms) {
      likelihood(private$.likelihood, private$.data, private$simulate, 
                 pars, initial.values, parms)
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
    #' @details 
    #' A mapping is a named argument, where name is the
    #' data colummn name, and value corresponds to the model variables (or an 
    #' expression to calculate from the model variables.)
    initialize = function(model, time, data, likelihood, ..., cumulative=FALSE, mapping=character()) {
      library(bbmle)
      if (!"Likelihood" %in% class(likelihood)) 
        stop("likelihood must be a Likelihood object")
      private$.likelihood = likelihood
      super$initialize(model, time, data, ..., cumulative = cumulative, mapping = mapping)
    }
  ),
  
  active = list(
    #' @field the names of all parameters
    parameters = function() {
      c(private$.model$parameters, private$.likelihood$par)
    }
  )
)
