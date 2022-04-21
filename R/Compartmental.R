#' An R6 class representing a compartmental model
#' 
#' This is a subclass of Model. Unlike a Model object which gives the
#' ODE for the rate of change of the states, a Compartmental model 
#' describes the model using compartments and transitions between 
#' compartments. Each transition may correspond to an event. The ODE 
#' system can then be derived by summing the transitions that flow from
#' the compartment (negative flows) and those that flow to the compartment
#' (positive flows). Because a compartmental model can describe events,
#' they can also be used for stochastic simulations (such as the Gillespie
#' method).
#' 
#' @docType class
#' @export
Compartmental <- R6Class(
  "Compartmental",
  inherit = Model,
  private = list(
    # the transitionss
    .transitions = list(),

    # add two rates
    add = function(a, b) {
      if (is.numeric(a)) {
        if (a == 0) return(b)
        if (is.numeric(b)) return(a+b)
      }
      if (is.numeric(b) && b == 0) a else
        as.call(list(as.name("+"), a, b))
    },
    
    # subtract two rates
    sub = function(a, b) {
      if (is.call(b)) {
        if (b[[1]] == "-" && length(b) == 2)
          return (private$add(a, b[[2]]))
      }
      if (is.numeric(a)) {
        if (is.numeric(b)) return(a-b)
        if (a == 0) return(as.call(list(as.name("-"), b)))
      }
      if (is.numeric(b) && b==0) return(a)
      as.call(list(as.name("-"), a, b))
    },
    
    #multiple a rate by a factor
    mul = function(a, b) {
      if (is.numeric(a)) {
        if (a == 0) return(0)
        if (a == 1) return(b)
        if (a == -1) return(private$sub(0, b))
        if (is.numeric(b)) return (a*b)
      }
      if (is.numeric(b)) {
        if (b == 0) return(0)
        if (b == 1) return(a)
        if (b == -1) return(private$sub(0, a))
      }
      if ((is.call(a) && a[[1]] == "/") || (is.call(b) && b[[1]] == "/")) {
        if (is.call(a) && a[[1]] == "/") {
          na = a[[2]]
          da = a[[3]]
        } else {
          na = a
          da = 1
        }
        if (is.call(b) && b[[1]] == "/") {
          nb = b[[2]]
          db = b[[3]]
        } else {
          nb = b
          db = 1
        }
        private$div(private$mul(na, nb), private$mul(da, db)) 
      } else if (is.call(a) && a[[1]] == "-" && length(a) == 2) {
        private$sub(0, private$mul(a[[2]], b))
      } else if (is.call(b) && b[[1]] == "-" && length(b) == 2) {
        private$sub(0, private$mul(a, b[[2]]))
      } else as.call(list(as.name("*"), a, b))
    },
    
    # divide a rate by a factor
    div = function(a, b) {
      if (is.numeric(b)) {
        if (b == 1) return(a)
        if (b == -1) return(private$sub(0, a))
        if (b == 0) stop("divide by 0")
      }
      if (is.call(a) && a[[1]] == "-" && length(a) == 2) {
        private$sub(0, private$div(a[[2]], b))
      } else if (is.call(b) && b[[1]] == "-" && length(b) == 2) {
        private$sub(0, private$div(a, b[[2]]))
      } else as.call(list(as.name("/"), a, b)) 
    },
    
    # formulate the equation for a compartment by adding the rates of
    # the transitions that flow to the compartment, and subtracting the
    # rates of the transitions that flow from the compartment.
    equation = function(compartment) {
      rate = Reduce(function(rate, y) {
        r = private$.formula[[y$name]]$formula
        if (!is.null(y$from) && y$from == compartment) {
          private$sub(rate, r) 
        } else if (!is.null(y$to) && y$to == compartment) {
          private$add(rate, r)
        } else rate
      }, private$.transitions, 0)
      super$compartment(call("~", as.name(compartment), rate))
    },

    # parse the formula to get the rate of a transition
    parse.rate = function(e) {
      if (!is.call(e)) as.character(e) else {
        if (e[[1]] != "~") stop("invalid transition") 
        list(compartment = e[[2]], rate=e[[3]])
      }  
    },
    
    # generate a name an unnamed transition
    transition.name = function(from, to) {
      name = paste0(from, "->", to)
      n = name
      i = 1
      while (!is.null(private$.transitions[[name]])) {
        name = paste0(n, ".", i)
        i = i + 1
      }
      name
    },
    
    # whether the name is defined as a transition, a compartment,
    # a parameter or a substitution
    defined = function(name) {
      super$defined || !is.null(private$.transitions[[name]])
    },
    
    # perform the actual rename. If from is a transition, the
    # rename the transition, otherwise, call the method from Model.
    do.rename = function(from, to) {
      if (!is.null(private$.transitions[[from]])) {
        private$.transitions[[to]] = private$.transitions[[from]]
        private$.transitions[[from]] = NULL
      } else super$do.rename(from, to)
    }
  ),
  
  public = list(
    #' @description constructor
    #' @param ... the names of the compartment, or substitutions
    #' @param title the name of the model
    #' @param file if not NULL, a path or connection to a model file 
    #' to read the model from
    #' @examples
    #' # an SIR model
    #' SIR = Compartmental$new(S, I, R, title="SIR")
    #' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
    #' SIR$transition(I->R ~ gamma*I, name="recovery")
    #' print(SIR)
    initialize = function(..., title = "", file = NULL) {
      self$title = title
      args = as.list(substitute(list(...)))[-1]
      ns = names(args)
      if (length(args) > 0) {
        for (i in 1:length(args)) {
          if (!is.null(ns) && ns[i] != "") {
            self$where(pairs=args[i])
          } else self$compartment(args[[i]])
        }
      }
    },
    
    #' @description define a compartment 
    #' @param name the name of the compartment
    #' @return an invisible self for chaining methods
    #' @examples 
    #' # an SIR model
    #' SIR = Compartmental$new(title="SIR")
    #' SIR$compartment("S")$compartment("I")$compartment(R)
    #' SIR$transition(S->I ~ beta*S*I/N, name="infection")
    #' SIR$transition(I->R ~ gamma*I, name="recovery")
    #' SIR$where(N=S+I+R)
    #' print(SIR)
    compartment = function(name) {
      if (length(name) == 0) {
        return(self)
      }
      if (!is.name(name) && !is.character(name))
        stop("invalid compartment name: ", name)
      else if (!is.name(name)) name = as.name(name)
      super$compartment(call("~", name, 0))
      invisible(self)
    },

    #' @description delete a compartment, a substitution or a transition
    #' @param name the name of the compartment, substitution or transition
    #' to delete.
    #' @return an invisible self for chaining methods
    #' @details The deleted compartment or substitution is converted to a 
    #' parameter. A parameter cannot be deleted. If the equation is changed
    #' so that a parameter is not used by the model any more, then the paramter
    #' is automatically removed.
    #' @examples
    #' # an SIR model
    #' SIR = Compartmental$new(S, I, R, title="SIR")
    #' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
    #' SIR$transition(I->R ~ gamma*I, name="recovery")
    #' print(SIR)
    delete = function(name) {
      if (!is.null(private$.transitions[[name]])) {
        self$transition(name=name)
      } else {
        if (!is.null(private$.compartments[[name]])) {
          for (tr in private$.transitions) {
            if (identical(tr$from, name) || identical(tr$to, name))
              self$transition(name=tr$name)
          }
        }
        super$delete(name)
      }
      invisible(self)
    },

    #' @description define or change a transition
    #'
    #' @param formula the transition can be fullt specified by a formula. Please see the
    #' details.
    #' @param from the name of the compartment that the transition originates
    #' @param to the name of the compartment that the transition flows into
    #' @param rate the rate of the transition, an expression. If rate is NULL, the transition is deleted.
    #' @param percapita whether the rate is per capita, i.e., the total rate would 
    #' be the per capita rate multiplied by the from compartment.
    #' @param name the name of the transition, a character. If NULL, the name
    #' @param ... named arguments specifying substitutions for the parameters
    #' used by this transition.
    #' is automatically generated
    #' @return the name of the transition
    #' @details The formula can be specified by a formula, which has the form
    #' from -> to ~ rate or to <- from ~ rate
    #' where the rate can be either an expression, or percapita(expression). 
    #' Alternatively the from, to, or rate (and percapita) can be specified by the
    #' arguments of the same name. 
    #' 
    #' A transition must be named. If the name is not provided, one is 
    #' automatically generated. Any further change to the transition must be
    #' referred by the name.
    #' 
    #' The from or to can be null. If from is NULL, then the transition is
    #' an input. If to is NULL, then the transition is an output. If both 
    # from and to are NULL, or rate is NULL, then the transition specified by
    #' the name is deleted.
    #' 
    #' If the transition with the given name exists, then the transition is changed.
    #' 
    #' @examples
    #' # an SIR model
    #' SIR = Compartmental$new(S, I, R, title="SIR")
    #' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
    #' SIR$transition(I->R ~ gamma*I, name="recovery")
    #' SIR$transition(S->R ~ v, percapita=TRUE, name="vaccination")
    #' # the following command changes the "vaccination" transition
    #' SIR$transition(S->NULL ~ v*S, name="vaccination")
    #' # this following command delete the "vaccination" transition, because
    #' # the rate is default to NULL.
    #' SIR$transition(name="vaccination")
    #' SIR$transition(NULL->S ~ lambda, name="births")
    #' SIR$transition(S -> NULL ~ percapita(mu), name="death.S")
    #' SIR$transition(I -> NULL ~ percapita(mu), name="death.I")
    #' SIR$transition(R -> NULL ~ percapita(mu), name="death.R")
    #' print(SIR)
    transition = function(formula = NULL, ..., from = NULL, to = NULL, rate = NULL, percapita = FALSE, name=NULL) {
      formula = substitute(formula)
      # parse the formula
      if (!is.null(formula)) {
        if (!is.call(formula) || formula[[1]] != "<-")
          stop("invalid transition")
        to = private$parse.rate(formula[[2]])
        if (is.null(to)) 
          stop("invalid transition")
        if (is.list(to)) {
          r = to$rate
          to = to$compartment
        } else r = NULL
        from = private$parse.rate(formula[[3]])
        if (is.list(from)) {
          if (!is.null(r)) 
            stop("invalid transition")
          r = from$rate
          from = from$compartment
        }
        if (!is.null(r)) rate = r
      } else rate = substitute(rate)

      # argument validity check
      if (!is.null(from) && !is.character(from)) {
        if (!is.name(from)) stop("invalid compartment ", from)
        from = as.character(from)
      }
      if (!is.null(to) && !is.character(to)) {
        if (!is.name(to)) stop("invalid compartment ", to)
        to = as.character(to)
      }
      if (!is.null(from) && is.null(private$.compartments[[from]]))
        stop("the compartment ", from, " is not defined")
      if (!is.null(to) && is.null(private$.compartments[[to]]))
        stop("the compartment ", to, " is not defined")
      
      # rate
      if (is.call(rate)) {
        if (rate[[1]] == "percapita") {
          rate = rate[[2]]
          percapita = TRUE
        }
      }
      if (percapita) 
        rate = private$mul(rate, as.name(from))
      # name
      if (is.null(name)) name = private$transition.name(from, to)

      # the transition exist?
      tr = private$.transitions[[name]]
      if (!is.null(tr)) {
        # delete the transition?
        private$.transitions[[name]] = NULL
        if (is.null(rate) || (is.null(from) && is.null(to))) {
          private$define.formula(name, NULL)
          private$.transitions[[name]] = NULL
          name = NULL
        } else { # change the transition
          rename = grepl(paste0(from, "->", to), name) &&
            (!identical(tr$from, from) || !identical(tr$to, to))
          new.name = if (rename) {
              private$transition.name(from, to)
          } else name
          private$define.formula(name, NULL)
          private$define.formula(new.name, rate)
          private$.transitions[[name]] = list(from = from, to = to, name=new.name)
          change.to = identical(to, tr$to)
          if (!identical(from, tr$from)) {
            if (!is.null(tr$from)) private$equation(tr$from)
          } else { #changed.to
            if (!is.null(tr$to)) private$equation(tr$to)
          }
          if (!is.null(from)) private$equation(from)
          if (!is.null(to)) private$equation(to)
          name = new.name
        }
      } else {
        private$define.formula(name, rate)
        private$.transitions[[name]] = list(from = from, to = to, name=name)
        if (!is.null(from)) private$equation(from)
        if (!is.null(to)) private$equation(to)
      }
      where = as.list(substitute(list(...)))[-1]
      if (length(where) > 0) for (i in 1:length(where)) {
        w = where[i]
        n = names(w)
        if (is.null(n) || n== "")
          stop("meaningless definition ", w)
        if (!is.null(private$.where[[n]]))
          stop("redefinition of ", n)
        self$where(pairs=w)
      }
      invisible(name)
    },
    
    #' @description format the class for printing
    format = function() {
      l = c(paste0("Compartmental: ", self$title))
      if (length(private$.compartments) > 0) {
        l = c(l, paste0("  Compartments: ", paste(self$compartments, collapse = ", ")))
      }
      if (length(private$.transitions) > 0) {
        l = c(l, "  transitions:")
        for (tr in private$.transitions)
          l = c(l, paste0("    \"", tr$name, "\" : ", 
                          if (is.null(tr$from)) "NULL" else tr$from, " -> ", 
                          if (is.null(tr$to)) "NULL" else tr$to,
                          " ~ ", deparse(private$.formula[[tr$name]]$formula)))
      }
      if (length(private$.where) > 0) {
        l = c(l, "  where")
        s = self$substitutions
        for (w in names(s))
          l = c(l, paste0("    ", w, " = ", deparse(s[[w]])))
      }
      if (length(self$parameters) > 0)
        l = c(l, paste0("  Parameters: ", paste(self$parameters, collapse=", ")))
      paste(l, collapse="\n")
    }
  ),
  
  active = list(
    #' @field transitions a read-only field to access the transitions.
    transitions = function() {
      lapply(
        private$.transitions,
        function(x) { 
          x$rate=private$.formula[[x$name]]$formula
          x
        }
      )
    },
    
    #' @field representation a read-only active field that returns the representation of the model
    #' it returns a list that contains self$compartments, self$transitions and
    #' self$substitutions. The compartmental model can then be reconstructed from the representation. 
    representation = function() {
      list(
        compartments = self$compartments,
        transitions = self$transitions,
        substitutions = self$substitutions
      )
    }
  )
)

if (exists("TEST") && is.logical(TEST) && TEST) {
  m = Compartmental$new(title="SIR", S, I, R)
  m$transition(I<-S ~ beta*S*I/N, N=S+I+R)
  m$transition(R<-I ~ gamma, percapita = TRUE)
  m$where(N=S+I+R)
  print(m)
}

