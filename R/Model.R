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
#' are called compartments, and the equations give the rate of change (i.e.,
#' the time derivatives of the states). An equation is specified by an
#' R formula, the parameters are automatically extracted. A model can then
#' be used to construct numerical or stochastic simulations, or be typeset
#' as latex equations.
#' 
#' @docType class
#' @export
Model <- R6Class(
  "Model",
  private = list(
    .compartments = list(),
    .where = list(),
    .formula = list(),
    .t = NULL,
    # external functions used by this model
    .external.functions = list(),
    
    # this function defines a formula with validity checking, 
    # and extract parameters is formula is NULL, then it is removed
    # if the formula with the given name already exists, it is redefined
    define.formula = function(name, formula) {
      # if formula is NULL, remove the formula
      if (is.null(formula)) {
        private$.formula[[name]] = NULL
        return()
      }
      e = if (is(formula, "Expression")) formula else Expression$new(formula)
      for (f in e$functions) {
        if (!is.null(private$.compartments[[f]]))
          stop("Redefining the compartment ", f, " as a function")
        if (!is.null(private$.where[[f]]))
          stop("Redefining the substitution ", f, " as a function")
        # check for parameter
        for (ff in private$.formula) {
          if (f %in% ff$parms)
            stop("Redefining the parameter ", ff, " as a function")
        }
      }
      private$.external.functions = union(private$.external.functions, e$functions)
      # are we redefining the formula or removing the formula?
      # if so, we remove the old definition
      private$.formula[[name]] = e
    },
    
    # perform a rename. The sanity check is already done in the rename method.
    do.rename = function(from, to) {
      update.defined = NULL
      if (!is.null(private$.compartments[[from]])) {
        formula.from = paste0(".d.", from)
        formula.to = paste0(".d.", to)
        private$.formula[[formula.to]] = private$.formula[[formula.from]]
        private$.formula[[formula.from]] = NULL
        private$.compartments[[to]] = private$.compartments[[from]]
        private$.compartments[[to]]$name = to
        private$.compartments[[to]]$value = formula.to
        private$.compartments[[from]] = NULL
        update.defined = TRUE
      } else if (!is.null(private$.where[[from]])) {
        formula.from = from
        formula.to = to
        private$.where[[to]] = private$.where[[from]]
        private$.where[[to]]$name = to
        private$.where[[to]]$value = to
        private$.where[[from]] = NULL
        change.formula = TRUE
        update.defined = TRUE
      } else stop(from, " is not defined")

      # change the name in each formula
      for (e in private$.formula) {
        e$rename(from, to)
      }
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
      private$.t = representation$.t
      if (!is.character(private$.t) || private$.t == "") private$.t = "t"
      for (C in names(representation$compartments)) {
        if (C == ".t") {
          private$.t = as.name(representation$compartments[[C]])
        } else {
          r = representation$compartments[[C]]
          self$compartment(call("~", as.name(C), r))
        }
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
    #' @param t the name of the independent variable, either a name or a string
    #' @param functions a list of functions to be used in this model
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
    initialize = function(..., t="t", functions=NULL, file=NULL, .restricted=FALSE) {
      self$restricted = .restricted
      if (!is.null(file)) private$load(file)
      private$.t = if (is.null(t) || t == "") "t" else if (is.character(t)) t else
        stop("Invalid independent variable name ", t)
      args = as.list(substitute(list(...)))[-1]
      if (!is.null(functions)) {
        nf = names(functions)
        if (is.null(nf) || any(nf == ""))
          stop("functions must be named")
      }
      if (length(functions) > 0) {
        for (i in 1:length(functions))
          self$attached.functions[[nf[i]]] = functions[[i]]
      }
      ns = names(args)
      if (length(args) > 0) {
        for (i in 1:length(args)) {
          if (!is.null(ns) && ns[i] != "") {
            self$where(pairs = args[i])
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
      if (name == private$.t)
        stop(name, " is the independent variable, and so cannot be used as a dependent variable name")
      if (!is.null(private$.where[[name]]))
        stop(name, " is already defined as a substitution, so cannot be used as a dependent variable name")
      formula = eq[[3]]
      if (is.null(formula)) {
        self$delete(name)
      } else {
        rate = paste0(".d.", name)
        private$define.formula(rate, formula)
        private$.compartments[[name]] = list(name = name, value = rate)
      }
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
        if (name == private$.t)
          stop(name, " is the independent variable, and so cannot be used as a parameter name")
        private$define.formula(name, defs[[name]])
        private$.where[[name]] = list(name=name, value=name)
      }
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
      if (name == private$.t)
        stop(name, " is the independent variable, and so cannot be deleted")
      if (!is.null(private$.compartments[[name]])) {
        private$.compartments[[name]] = NULL
        private$define.formula(paste0(".d.", name), NULL)
      } else if (!is.null(private$.where[[name]])) {
        private$.where[[name]] = NULL
        private$define.formula(name, NULL)
      } else stop(name, " is not defined")
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
    #' # an SIR model
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
      if (from == private$.t) {
        private$.t = as.name(from)
        return(invisible(self))
      }
      
      if (is.null(to))
        return(self$delete(from))
      if (!is.name(to) && ! is.character(to))
        stop("Invalid name ", to)
      to = as.character(to)
      if (to == private$.t)
        stop(from, " cannot be renamed to the independent variable")

      if (from == to) return()

      private$do.rename(from, to)
      invisible(self)
    },

    #' @description format the class for printing
    format = function() {
      l = c("Model:", paste("independent variable:", private$.t))
      if (length(private$.compartments) > 0) {
        l = c(l, "  Compartments:")
        for (c in private$.compartments) {
          f = private$.formula[[c$value]]
          l = c(l, paste0("    ", c$name, " ~ ", deparse(f$expr)))
        }
      }
      if (length(private$.where) > 0) {
        l = c(l, "  where")
        for (w in private$.where) {
          l = c(l, paste0("    ", w$name, " = ", deparse(private$.formula[[w$value]]$expr)))
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
          call("==", call("'", as.name(C$name)), private$.formula[[C$value]]$expr)
        }
      )
      where = lapply(
        private$.where,
        function(w) {
          call("==", as.name(w$name), private$.formula[[w$value]]$expr)
        }
      )
      list(equations = compartments, where = where)
    },
    
    #' @field parameters A read-only field that returns a character vector of parameter names
    parameters = function() {
      parms = c()
      for (e in private$.formula) {
        if (length(e$parms) > 0) {
          parms = union(parms, e$parms)
        }
      }
      remove = c(names(private$.compartments), names(private$.where), private$.t)
      setdiff(parms, remove)
    },

    #' @field substitutions A read-only field that returns a named list of expressions
    substitutions = function() {
      l = list()
      compartments = self$compartments
      for (w in private$.where) {
        e = private$.formula[[w$value]]
        l[[w$name]] = e$expr
        if (length(intersect(e$parms, compartments)) > 0)
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
        .t = private$.t,
        compartments = sapply(
          private$.compartments,
          function(C) private$.formula[[C$value]]$expr
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
    },
    
    #' @field t the independent variable anme
    t = function() {
      private$.t
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