# ==============================================================================
# RGillespie: Gillespie SSA (pure R implementation)
# ==============================================================================
#
# This file implements a Simulator subclass that runs a Gillespie / SSA style
# stochastic simulation for a Compartmental model.
#
# Design:
#   - The model is compiled into a function:
#       f(time, y, parms) -> c(new_time, new_y, extra_outputs)
#   - One call to f() advances the system by exactly one event:
#       1) compute all transition rates
#       2) sample the next event type proportional to its rate
#       3) sample the waiting time ~ Exp(total_rate)
#       4) update y according to the selected transition (±1 in compartments)
#
# Notes:
#   - This implementation assumes *integer* compartment counts.
#   - Transition rates are evaluated using the same substitution mechanism as
#     deterministic simulators (bind y/parms, evaluate substitutions).
#   - Substitutions are evaluated every event (simple and safe, but can be
#     optimized later if needed).
# ==============================================================================

#' R6 class implementing the Gillespie method in R
#'
#' This is a subclass of [Simulator], using a pure R implementation of the
#' Gillespie stochastic simulation algorithm to simulate a [Compartmental] model.
#'
#' @docType class
#' @examples
#' # an SIR model
#' SIR <- Compartmental$new(S, I, R)
#' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
#' SIR$transition(I->R ~ gamma*I, name="recovery")
#' g <- RGillespie$new(SIR)
#' g$simulate(0:100, y0=c(S=1000, I=10, R=0), parms=c(beta=0.4, gamma=0.2))
#' @export
RGillespie <- R6Class(
  "RGillespie",
  inherit = Simulator,
  
  private = list(
    # --------------------------------------------------------------------------
    # Small code-generation helpers
    # --------------------------------------------------------------------------
    
    # Format the transition rate vector as an R call: c(rate1, rate2, ...)
    format.rate = function(T) {
      as.call(c(quote(c), lapply(unname(T), function(tr) tr$rate)))
    },
    
    # Format substitutions as assignment calls: x <- expr
    format.substitution = function(S) {
      lapply(names(S), function(var) call("<-", as.name(var), S[[var]]))
    },
    
    # Bind named entries of y/parms into variables of the same names.
    #
    # Produces:  S <- y[[1]] ; I <- y[[2]] ; ...
    format.var = function(S, name) {
      lapply(S, function(var) {
        call("<-", as.name(var), call("[[", as.name(name), match(var, S)))
      })
    },
    
    # --------------------------------------------------------------------------
    # Compile model -> single-event step function(time, y, parms)
    # --------------------------------------------------------------------------
    build = function(model) {
      alias_names <- names(private$alias)
      extra_syms <- lapply(alias_names, as.name)
      names(extra_syms) <- alias_names
      
      # Compiled function body:
      # {
      #   <bind y/parms>
      #   <compute substitutions>
      #   .rates <- c(rate_1, ..., rate_M)
      #   .sum   <- cumsum(.rates)
      #   .total <- .sum[[M]]
      #   if (!is.finite(.total) || .total == 0) return(c(Inf, y))
      #   .transition <- which(.sum >= runif(1) * .total)[[1]]
      #   time <- time + rexp(1, .total)
      #   switch(.transition, { update y }, { update y }, ...)
      #   c(time, y, <extras>)
      # }
      body_calls <- c(
        as.name("{"),
        private$format.var(model$compartments, "y"),
        private$format.var(model$parameters, "parms"),
        private$format.substitution(model$substitutions),
        
        call("<-", quote(.rates), private$format.rate(model$transitions)),
        quote(.sum <- cumsum(.rates)),
        call("<-", quote(.total), call("[[", quote(.sum), length(model$transitions))),
        quote(if (!is.finite(.total) || .total == 0) return(c(Inf, y))),
        
        quote(.transition <- which(.sum >= runif(1) * .total)[[1]]),
        call("<-", as.name(model$t), call("+", as.name(model$t), quote(rexp(1, .total)))),
        
        # switch over transitions: each case updates y by ±1
        as.call(c(
          quote(switch),
          quote(.transition),
          unname(lapply(model$transitions, function(tr) {
            updates <- list()
            
            if (!is.null(tr$from)) {
              idx_from <- which(tr$from == model$compartments)
              updates <- c(updates, call(
                "<-",
                call("[[", quote(y), idx_from),
                call("-", as.name(tr$from), 1)
              ))
            }
            
            if (!is.null(tr$to)) {
              idx_to <- which(tr$to == model$compartments)
              updates <- c(updates, call(
                "<-",
                call("[[", quote(y), idx_to),
                call("+", as.name(tr$to), 1)
              ))
            }
            
            as.call(c(as.name("{"), updates))
          }))
        )),
        
        as.call(c(list(as.name("c"), as.name(model$t), as.name("y")), extra_syms))
      )
      
      args <- alist(,,)
      names(args) <- c(model$t, "y", "parms")
      as.function(c(args, as.call(body_calls)))
    },
    
    # --------------------------------------------------------------------------
    # Run the Gillespie simulation on a time grid t
    # --------------------------------------------------------------------------
    .simulate = function(t, y, parms) {
      # Validate integer, nonnegative initial conditions
      mistype <- which(is.na(y) | y != as.integer(y) | y < 0)
      if (length(mistype) != 0) {
        s <- if (length(mistype) > 1) "s" else ""
        a <- if (length(mistype) > 1) "" else "a "
        stop(
          "the initial condition", s, " for ",
          paste(names(y[mistype]), collapse = ", "),
          " must be ", a, "nonnegative integer", s, "."
        )
      }
      
      y <- y[private$compartments]
      parms <- parms[private$parameters]
      
      if (any(!is.finite(t)))
        stop("invalid time")
      if (length(t) < 2)
        stop("need at least two time points")
      if (is.null(names(y)) || any(names(y) == ""))
        stop("the initial values must be named")
      
      # Output storage
      data <- data.frame(
        time = t,
        matrix(NA_real_, nrow = length(t), ncol = length(y))
      )
      colnames(data) <- c("time", names(y))
      data[1, -1] <- y
      
      i <- 2
      time <- t[1]
      
      while (TRUE) {
        step <- self$model(time, y, parms)
        time <- step[[1]]
        
        # Fill output rows until the next event time crosses them
        while (time > t[[i]]) {
          data[i, -1] <- y
          i <- i + 1
          if (i > length(t)) {
            return(data)
          }
        }
        
        # Update state vector (step returns c(time, y, extras...); we ignore extras here)
        y <- step[2:(1 + length(y))]
      }
    }
  )
)
