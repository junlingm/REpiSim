#' R6 class implementing the Gillespie method in R
#'  
#' This is a subclass of Simulator, using an R implementation of the Gillespie method
#' to simulate a compartmental model.
#' 
#' @docType class
#' @examples
#' # an SIR model
#' SIR = Compartmental$new(S, I, R, title="SIR")
#' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
#' SIR$transition(I->R ~ gamma*I, name="recovery")
#' g = RGillespie$new(SIR)
#' g$simulate(0:100, y0=c(S=1000, I=10, R=0), parms=c(beta=0.4,gamma=0.2))
#' @export
RGillespie = R6Class(
  "RGillespie",
  inherit = Simulator,
  private = list(
    # format the transition rates as R commands
    format.rate = function(T) {
      as.call(c(quote(c), lapply(unname(T), function(tr) tr$rate)))
    },

    # format the substitutions as R commands
    format.substitution = function(S) {
      lapply(names(S), function(var) {
        call("<-", as.name(var), S[[var]])
      })
    },
    
    # format assign a compenent of a named vector as a variable of the same name 
    format.var = function(S, name) {
      lapply(S, function(var) {
        call("<-", as.name(var), call("[[", as.name(name), which(var == S)))
      })
    },
    
    # return the R code for simulating the model for a single time step (one event)
    build = function(model) {
      l = c(
        as.name("{"),
        private$format.var(model$compartments, "y"),
        private$format.var(model$parameters, "parms"),
        private$format.substitution(model$substitutions),
        call(
          "<-", quote(.rates), 
          private$format.rate(model$transitions)
        ),
        quote(.sum <- cumsum(.rates)),
        call("<-", quote(.total), call("[[", quote(.sum), length(model$transitions))),
        quote(if (!is.finite(.total) || .total == 0) return(c(Inf, y))),
        quote(.transition <- which(.sum >= runif(1) * .total)[[1]]),
        quote(t <- t + rexp(1, .total)),
        as.call(c(
          quote(switch),
          quote(.transition),
          unname(lapply(model$transitions, function(tr) {
            l = list()
            if (!is.null(tr$from))
              l = c(l, call(
                "<-",
                call("[[", quote(y), which(tr$from == model$compartments)),
                call("-", as.name(tr$from), 1)
              ))
            if (!is.null(tr$to))
              l = c(l, call(
                "<-",
                call("[[", quote(y), which(tr$to == model$compartments)),
                call("+", as.name(tr$to), 1)
              ))
            as.call(c(as.name("{"), l))
          })
        ))),
        quote(c(t, y))
      )
      as.function(c(alist(t=, y=, parms=), as.call(l)))
    },
    
    # run the simulation
    run = function(t, y, parms) {
      y = y[private$compartments]
      parms = parms[private$parameters]      
      if (any(!is.finite(t)))
        stop("invalid time")
      if (length(t) < 2)
        stop("need at least two time points")
      if (is.null(names(y)) || any(names(y) == ""))
        stop("the initial values must be named")
      data = data.frame(
        t,
        matrix(NA, nrow=length(t), ncol=length(y))
      )
      colnames(data)= c("time", names(y))
      data[1, -1] = y
      i = 2
      time = t[1]
      while (TRUE) {
        l = self$model(time, y, parms)
        time = l[[1]]
        while (time > t[[i]]) {
          data[i, -1] = y
          i = i + 1
          if (i > length(t)) {
            return(data)
          }
        }
        y = l[-1]
      }
    }
  )
)
