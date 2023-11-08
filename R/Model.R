D = function(x, ...) {
  vars = as.list(substitute(list(...)))[-1]
  as.call(c(list("'", x), vars))
}

letter.label = function(n, enclosed = c("(", ")")) {
  if (!is.null(enclosed) && length(enclosed) == 1)
    enclosed = rep(enclosed, 2)
  paste0(enclosed[[1]], letters[n], enclosed[[2]])
}

#' R6 class representing a mathematical model.
#' 
#' A mathematical model is described by a system of ODEs. The state variables
#' are called compartments, and the equations gives the rate of change (i.e.,
#' the time derivatives of the states). The equations are specified by
#' R formula, the parameters are automatically extracted. A model can then
#' be used to construct numerical simulations or stochastic simulations,
#' or be formatted as latex equations.
#' 
#' @docType class
#' @export
Model <- R6Class(
  "Model",
  private = list(
    .compartments = list(),
    .parameters = list(),
    .where = list(),
    .formula = list(),
    # whether the compartment field of a substitution need to be recalculated
    recalc.compartment = FALSE,
    # external functions used by this model
    .external.functions = list(),

    # build up the order of dependence of an alias
    build.order = function(order, info) {
      # if var is already available in order, no need to change the order.
      if (info$name %in% order) return(order)
      # calculate the unavailable dependencies (i.e., not in order)
      deps = setdiff(private$.formula[[info$value]]$depend, order)
      # no dependencies or dependencies are already available
      # then we can just put it in order
      if (length(deps) == 0) return(c(order, info$name))
      for (d in deps) {
        v = private$.where[[d]]
        if (!is.null(v)) order = private$build.order(order, v)
      }
      c(order, info$name)
    },
    
    # this function returns an alist with a given number n of arguments
    make.alist = function(n) {
      if (n == 0) NULL else {
        l = alist(a=)
        for (i in 1:n) {
          names(l)[i] = letter.label(i, enclose=NULL)
        }
        c(l, NULL)
      }
    },

    # perform a substitute that replaces a given variable name 
    # in the formula by an expression, given in the named list subs
    substitute = function(formula, subs) {
      if (is.name(formula)) {
        name = as.character(formula)
        return(if (!is.null(subs[[name]])) subs[[name]] else formula)
      }
      if (is.call(formula)) {
        name = as.character(formula[[1]])
        if (!is.null(subs[[name]])) 
          formula[[1]] = as.name(subs[[name]])
        n = length(as.list(formula))
        for (i in 2:n)
          formula[[i]] = private$substitute(formula[[i]], subs)
      }
      formula
    },
    
    # this function extracts all parameters (names) used in the 
    # given formula which name is name. It returns a named list, where the name is 
    # the variable name, and the value is a list containing the 
    # name and if the name is a function, the alist of the function
    # definition. The argument parameters gives the initial value 
    # where an empty list means that no parameters have been found yet.
    extract.parameters = function(name, formula, parameters = list()) {
      new = function(var, definition = NULL) {
        p = private$.parameters[[var]]
        if (is.null(p)) {
          list(name=var, defined = name, definition = definition) 
        } else {
          p$defined = c(p$defined, name)
          p
        }
      }
      if (is.name(formula)) {
        var = as.character(formula)
        if (var != "t" && is.null(parameters[[formula]]))
          parameters[[var]] = new(var)
      } else if (is.call(formula)) {
        var = as.character(formula[[1]])
        if (! var %in% private$.external.functions)
          private$.external.functions = c(private$.external.functions, var)
        for (t in as.list(formula)[-1])
          parameters = private$extract.parameters(name, t, parameters)
      }
      parameters
    },
    
    # when a formula which name is given by owner is deleted, 
    # the parameters that is used in it will have the dependence 
    # on the formula removed. if there it is used in no formula
    # then the parameter will be deleted.
    remove.owner = function(l, owner) {
      Filter(
        function(x) !is.null(x),
        lapply(
          l,
          function(p) {
            if (owner %in% p$defined) {
              p$defined = p$defined[p$defined != owner]
              if (length(p$defined) == 0) NULL else p
            } else p
          }
        )
      )
    },
    
    # this function defines a formula with validity checking, 
    # and extract parameters is formula is NULL, then it is removed
    # if the formula with the given name already exists, it is redefined
    define.formula = function(name, formula) {
      pars = private$extract.parameters(name, formula)
      for (p in pars) {
        var = p$name
        if (!is.null(p$definition)) {
          if (!is.null(private$.compartments[[var]]))
            stop("Redefining the compartment ", var, " as a function")
          if (!is.null(private$.where[[var]]))
            stop("Redefining the substitution ", var, " as a function")
          if (!is.null(private$.parameters[[var]]))
            stop("Redefining the parameter ", var, " as a function")
        }
      }
      # are we redefining the formula or removing the formula?
      # if so, we remove the old definition
      if (!is.null(private$.formula[[name]]) || is.null(formula)) {
        private$.parameters = private$remove.owner(private$.parameters, name)
        if (is.null(formula)) return()
      }
      private$.formula[[name]] = list(formula = formula)
      # define parameters
      for (p in pars) {
        var = p$name
        if (!is.null(private$.compartments[[var]]) || !is.null(private$.where[[var]])) {
          if (var != name)
            private$.formula[[name]]$depend = c(private$.formula[[name]]$depend, var)
        } else private$.parameters[[var]] = p
      }
    },
    
    # whether the name is defined as a compartment, a substitution or a parameter
    defined = function(name) {
      !is.null(private$.compartments[[name]]) || 
      !is.null(private$.where[[name]]) ||
      !is.null(private$.parameters[[name]])
    },

    # perform a rename. The sanity check is already done in the rename method.
    do.rename = function(from, to) {
      if (!is.null(private$.compartments[[from]])) {
        formula.from = paste0(".d.", from)
        formula.to = paste0(".d.", to)
        private$.compartments[[to]] = private$.compartments[[from]]
        private$.compartments[[to]]$name = to
        private$.compartments[[to]]$value = formula.to
        private$.compartments[[from]] = NULL
        private$recalc.compartment = length(private$.where) > 0
      } else if (!is.null(private$.where[[from]])) {
        formula.from = from
        formula.to = to
        private$.where[[to]] = private$.where[[from]]
        private$.where[[to]]$name = to
        private$.where[[to]]$value = to
        private$.where[[from]] = NULL
        change.formula = TRUE
        private$recalc.compartment = length(private$.where) > 0
      } else if (!is.null(private$.parameters[[from]])) {
        formula.from = NULL
        private$.parameters[[to]] = private$.parameters[[from]]
        private$.parameters[[to]]$name = to
        private$.parameters[[from]] = NULL
      } else stop(from, " is not defined")
      if (!is.null(formula.from)) {
        private$.formula[[formula.to]] = private$.formula[[formula.from]]
        private$.formula[[formula.from]] = NULL
        # if to was a parameter, remove the parameter
        if (!is.null(private$.parameters[[to]])) 
          private$.parameters[[to]] = NULL
      }
      
      # change the name in each formula
      subs = list()
      subs[[from]] = as.name(to)
      private$.formula = lapply(
        private$.formula,
        function(f) {
          f$formula = private$substitute(f$formula, subs=subs)
          f$depend[f$depend == from] = to
          f
        }
      )
      # change the defined field of the parameters
      private$.parameters = lapply(
        private$.parameters,
        function(p) {
          p$defined[p$defined == from] = to
          p
        }
      )
    },
    
    # calculate if the substitution s depends on a compartment
    calc.compartment = function(s) {
      f = private$.formula[[s$value]]
      for (d in f$depend) {
        if (!is.null(private$.compartments[[d]])) 
          return (TRUE)
        w = private$.where[[d]]
        if (!is.null(w)) {
          if (is.null(w$compartment))
            private$.where[[d]]$compartment = private$calc.compartment(w)
          if (private$.where[[d]]$compartment) return (TRUE)
        }
      }
      s$compartment = FALSE
      s
    },
    
    # load the model from a file
    load = function(file) {
      if (is.character(file)) {
        if (!file.exists(file))
          stop("file does not exist: ", file)
        e = new.env()
        load(file, envir=e)
        file = e$model
        if (is.null(file)) 
          stop("not a valid model file: ", file)
      }
      private$construct(file)
      self$attached.functions = file$attached.functions
    },
    
    #reconstruct the model from a representation
    construct = function(representation) {
      if (!identical(representation$class, "Model"))
        stop("invalid model file")
      for (C in names(representation$compartments)) {
        r = representation$compartments[[C]]
        self$compartment(call("~", as.name(C), r))
      }
      self$where(pairs=representation$substitutions)
    }
  ),

  public = list(
    #' @field restricted a boolean variable indicating whether to restrict functions 
    #' allowed to be used in the ODE system. Default to FALSE
    restricted = FALSE,
    
    #' @field attached.functions the list of provided R functions to be used in
    #' the model
    attached.functions = list(),
    
    #' @description constructor
    #' 
    #' It constructs a Model object with compartments and substitutions.
    #' 
    #' @param ... Each extra parameter is either passed to the `compartment` method
    #' if it is a formula with the form `name ~ value`, or passed to the
    #' `where` methods to define a substitution if it is a named expression.
    #' @param file if not NULL, a path or connection to a model file 
    #' to read the model from
    #' @param .restricted a boolean variable indicating whether the ODE system
    #' only has access to a limited set of functions. Default to TRUE
    #' @examples
    #' # An SIR model 
    #' SIR = Model$new(
    #'   S ~ -beta*S*I/N, # the dS/dt equation
    #'   I ~ beta*S*I/N - gamma*I, # the dI/dt equation
    #'   R ~ gamma*I, # the dR/dt equation
    #'   N = S + I + R # the total population N
    #' )
    #' print(SIR)
    initialize = function(..., file=NULL, .restricted=FALSE) {
      self$restricted = .restricted
      if (!is.null(file)) private$load(file)
      args = as.list(substitute(list(...)))[-1]
      ns = names(args)
      if (length(args) > 0) {
        for (i in 1:length(args)) {
          if (!is.null(ns) && ns[i] != "") {
            if (is.call(args[[i]]) && args[[i]][[1]] == "function") {
              self$attached.functions[[ns[i]]] = eval(args[[i]])
            } else self$where(pairs = args[i])
          } else self$compartment(args[[i]])
        }
      }
    },
    
    #' @description
    #' define a compartment using a formula
    #' 
    #' @param eq a formula  in the form of `name ~ rate`. This defined a 
    #' compartment with a given name and rate of change.
    #' @return The invisible object `self` for chained operations
    #' 
    #' @details The compartment name can coincide with a parameter name, 
    #' in which case the parameter is converted into a compartment. But the
    #' name cannot conflict with another compartment or a substitution.
    #' 
    #' @examples 
    #' # an SIR model
    #' SIR = Model$new()
    #' SIR$compartment(S ~ -beta*S*I/N)$
    #'   compartment(I ~ beta*S*I/N - gamma*I)$
    #'   compartment(R ~ gamma*I)$
    #'   where(N = S + I + R)
    #' print(SIR)
    compartment = function(eq) {
      if (!is.call(eq) || eq[[1]] != "~")
        stop("Invalid equation")
      name = as.character(eq[[2]])
      formula = eq[[3]]
      if (is.null(formula)) {
        self$delete(name)
      } else {
        rate = paste0(".d.", name)
        private$define.formula(rate, formula)
        if (is.null(private$.compartments[[name]]))
          private$.compartments[[name]] = list(name = name, value = rate)
        p = private$.parameters[[name]]
        if (!is.null(p)) {
          for (d in p$defined) {
            private$.formula[[d]]$depend = c(private$.formula[[d]]$depend, name)
          }
          private$.parameters[[name]] = NULL
        }
      }
      private$recalc.compartment = length(private$.where) > 0
      invisible(self)
    },
    
    #' @description
    #' define parameters as substitutions
    #' 
    #' @param ... Each extra parameter must be named, which name is a parameter
    #' and the value is the substitution.
    #' @param pairs a named list, which is an alternative form to provide substitutions
    #' 
    #' @details If the expression uses a compartment, then the parameter is considered
    #' a compartment when simulating the model, i.e., it is returned as a column
    #' in the data frame returned by the simulation.
    #' 
    #' @examples
    #' # an SIR model
    #' SIR = Model$new()
    #' SIR$compartment(S ~ -beta*S*I)$
    #'   compartment(I ~ beta*S*I - gamma*I)$
    #'   compartment(R ~ gamma*I)$
    #'   where(pairs=list(beta=quote(b/N)))$
    #'   where(N = S + I + R)
    #' print(SIR)
    where = function(..., pairs = NULL) {
      # combine .... and pairs
      defs = c(as.list(substitute(list(...)))[-1], pairs)
      # if called without an argument, return all substitutions
      if (length(defs) == 0) 
        return(self$substitutions)

      # sanity checks
      ns = names(defs)
      if (is.null(ns) || any(ns == "")) 
        stop("Substitutions must be in the form name=value")
      u = unique(ns)
      if (length(u) != length(ns))
        stop("Redefined substitution",
             if (length(u) - length(ns) > 1) "s " else " ",
             paste(unique(setdiff(ns, u)), collapse=", ")
        )

      compartments = self$compartments
      redef = ns[ns %in% compartments]
      if (length(redef) > 0)
        stop("Redefining compartment", 
             if (length(redef) > 1) "s " else " ",
             paste(redef, collapse = ", ")
        )

      # define each substitution
      for (name in ns) {
        private$define.formula(name, defs[[name]])
        if (is.null(private$.where[[name]]))
          private$.where[[name]] = list(name=name, value=name)
        p = private$.parameters[[name]]
        if (!is.null(p)) {
          for (d in p$defined) {
            private$.formula[[d]]$depend = c(private$.formula[[d]]$depend, name)
          }
          private$.parameters[[name]] = NULL
        }
      }
      private$recalc.compartment = length(private$.where) > 0
      invisible(self)
    },

    #' @description
    #' delete a compartment or a substitution
    #' 
    #' @param name a character specifying the  name of the compartment 
    #' or the substitution to be deleted.
    #' @return an invisible self to chain methods.
    #' 
    #' @details The deleted compartment or substitution is converted to a 
    #' parameter. A parameter cannot be deleted. If the equation is changed
    #' so that a parameter is not used by the model any more, then the paramter
    #' is automatically removed.
    #' 
    #' @examples
    #' # an SIR model
    #' SIR = Model$new()
    #' SIR$compartment(S ~ -beta*S*I/N)$
    #'   compartment(I ~ beta*S*I/N - gamma*I)$
    #'   compartment(R ~ gamma*I)$
    #'   where(N = S + I + R)
    #' SIR$delete("N")
    #' print(SIR)
    delete = function(name) {
      if (length(name) > 1) {
        for (n in name) self$delete(n)
        return(invisible(self))
      }
      if (!is.name(name) && !is.character(name))
        stop("invalid name ", name)
      if (!is.null(private$.compartments[[name]])) {
        private$.compartments[[name]] = NULL
        private$define.formula(paste0(".d.", name), NULL)
      } else if (!is.null(private$.where[[name]])) {
        private$.where[[name]] = NULL
        private$define.formula(name, NULL)
      } else stop(name, " is not defined")
      # if name is used, make it a parameter
      used = Filter(
        function(f) { name %in% f$depend },
        private$.formula
      )
      private$.parameters[[name]] = if (length(used) == 0) NULL else
        list(name = name, defined = used)
      # remove name from the dependences of each formula
      private$.formula = lapply(
        private$.formula,
        function(f) {
          f$depend = f$depend[f$depend != name]
          f
        }
      )
      private$recalc.compartment = length(private$.where) > 0
      invisible(self)
    },

    #' @description rename a compartment, a substitution or a parameter
    #'
    #' @param from is the name of the compartment, substitution or paramter
    #' to be changed.
    #' @param to is the new name
    #' @param formula the @param from or @param to can be specified by a formula
    #' in the form of `from -> to` or `to <- from`.
    #' @return an invisible self to chain methods
    #' 
    #' @details The new name cannot conflict with another compartment,
    #' substitution or paramter
    #' 
    #' @examples
    #' an SIR model
    #' SIR = Model$new()
    #' SIR$compartment(S ~ -beta*S*I/N)$
    #'   compartment(I ~ beta*S*I/N - gamma*I)$
    #'   compartment(R ~ gamma*I)$
    #'   where(N = S + I + R)
    #' SIR$rename(S->U)$rename(beta->b)
    #' print(SIR)
    rename = function(formula=NULL, from=NULL, to=NULL) {
      formula = substitute(formula)
      if (!is.null(formula)) {
        if (!is.call(formula) || formula[[1]] != "<-")
          stop("invalid formula")
        to = formula[[2]]
        from = formula[[3]]
      }

      if (is.null(from))
        stop("The name to be changed is missing")
      if (!is.name(from) && ! is.character(from))
        stop("Invalid name ", from)
      from = as.character(from)
      if (!private$defined(from))
        stop(from, " is not defined")

      if (is.null(to))
        return(self$delete(from))
      if (!is.name(to) && ! is.character(to))
        stop("Invalid name ", to)
      to = as.character(to)
      if (private$defined(to))
        stop(to, " is already defined")
      
      if (from == to) return()

      private$do.rename(from, to)
      invisible(self)
    },

    #' @description format the class for printing
    format = function() {
      l = "Model: "
      if (length(private$.compartments) > 0) {
        l = c(l, "  Compartments:")
        for (c in private$.compartments) {
          f = private$.formula[[c$value]]
          l = c(l, paste0("    ", c$name, " ~ ", deparse(f$formula)))
        }
      }
      if (length(private$.where) > 0) {
        l = c(l, "  where")
        for (w in private$.where) {
          l = c(l, paste0("    ", w$name, " = ", deparse(private$.formula[[w$value]]$formula)))
        }
      }
      if (length(self$parameters) > 0)
        l = c(l, paste0("  Parameters: ", paste(self$parameters, collapse=", ")))
      paste(l, collapse="\n")
    }
  ),
  
  active = list(
    #' @field compartments A read-only field that returns a character vector of compartment names
    compartments = function() {
      names(private$.compartments)
    },

    #' @field equations A read-only field that returns a named list model equations
    equations = function() {
      compartments = lapply(
        private$.compartments,
        function(C) {
          call("==", call("'", as.name(C$name)), private$.formula[[C$value]]$formula)
        }
      )
      where = lapply(
        private$.where,
        function(w) {
          call("==", as.name(w$name), private$.formula[[w$value]]$formula)
        }
      )
      list(equations = compartments, where = where)
    },
    
    #' @field parameters A read-only field that returns a character vector of parameter names
    parameters = function() {
      names(private$.parameters)
    },

    #' @field substitutions A read-only field that returns a named list of expressions
    substitutions = function() {
      l = list()
      compartments = self$compartments
      if (private$recalc.compartment) {
        private$.where = lapply(private$.where, function(w) {
          w$compartment = NULL
          w
        })
        private$.where = lapply(private$.where, function(w) {
          w$compartment = private$calc.compartment(w)
          w
        })
        private$recalc.compartment = FALSE
      }
      for (w in private$.where) {
        if (!is.null(w$compartment))
          l[[w$name]] = private$.formula[[w$value]]$formula
          attr(l[[w$name]], "compartment") = TRUE
      }
      l
    },
    
    #' @field representation a read-only active field that returns the representation of the model
    #' it returns a list that contains the equations substitutions. The compartmental model can 
    #' then be reconstructed from the representation. 
    representation = function() {
      list(
        class = "Model",
        compartments = sapply(
          private$.compartments,
          function(C) private$.formula[[C$value]]$formula
        ),
        substitutions = self$substitutions,
        attached.functions = self$attached.functions
      )
    },
    
    #' @field order a read-only active field that returns the order of the equations
    #' and aliases that must appear to satisfy dependencies in calculation
    order = function() {
      Reduce(private$build.order, private$.where, c())
    },
    
    #' @field missing a read-only field that returns the names of functions defined 
    #' neither in the global environment, nor in the attached.functions list.
    missing = function() {
      l = c()
      for (fname in private$.external.functions) {
        f = if (exists(fname)) get(fname) else self$attached.functions[[fname]]
        if (!is.function(f))
          l = c(l, fname)
      }
      l
    }
  )
)

if (exists("TEST") && is.logical(TEST) && TEST) {
  m = Model$new(
    S ~ - beta*S*I/N,
    I ~ beta*S*I/N - gamma*I,
    R ~ gamma*I,
    N = S+I+R
  )
  print(m)
}