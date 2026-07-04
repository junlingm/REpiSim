# ==============================================================================
# C++ expression and model-context converter
# ==============================================================================
#
# Internal helper used by compiled simulators. It converts the symbolic pieces
# shared by simulator backends into C++ fragments:
#   - R expressions -> C++ expressions
#   - state/parameter vector entries -> named local variables
#   - substitutions -> ordered local assignments
#   - attached R functions -> cpp11::function bindings
#
# Simulator-specific code should still own its own driver semantics, such as
# Gillespie event selection or ODE derivative output.
# ==============================================================================

CppConverter <- R6Class(
  "CppConverter",
  
  public = list(
    compartments = NULL,
    parameters = NULL,
    alias = NULL,
    
    initialize = function(compartments, parameters, alias) {
      self$compartments <- compartments
      self$parameters <- parameters
      self$alias <- alias
    },
    
    remove_bracket = function(x) {
      if (is.call(x) && identical(x[[1]], as.name("(")))
        self$remove_bracket(x[[2]])
      else x
    },
    
    expr = function(C) {
      if (is.call(C)) {
        # normalize R's [[ to [ for our purposes
        if (identical(C[[1]], as.name("[["))) C[[1]] <- as.name("[")
        
        op <- as.character(C[[1]])
        
        switch(
          op,
          `+` = if (length(C) == 2) self$expr(C[[2]]) else
            paste0(self$expr(C[[2]]), " + ", self$expr(C[[3]])),
          `-` = if (length(C) == 2) {
            paste0("-", self$expr(C[[2]]))
          } else {
            paste0(self$expr(C[[2]]), " - ", self$expr(C[[3]]))
          },
          `*` = paste0(self$expr(C[[2]]), " * ", self$expr(C[[3]])),
          `/` = paste0(self$expr(C[[2]]), " / ", self$expr(C[[3]])),
          `^` = paste0(
            "pow(",
            self$expr(self$remove_bracket(C[[2]])),
            ", ",
            self$expr(self$remove_bracket(C[[3]])),
            ")"
          ),
          `(` = paste0("(", self$expr(C[[2]]), ")"),
          `[` = if (length(C) > 3) {
            # support nested indexing like x[i][j] by recursively emitting calls
            self$expr(C[-1])
          } else {
            paste0(self$expr(C[[2]]), "[", self$expr(C[[3]]), "]")
          },
          {
            # generic function call
            paste0(
              "cpp11::as_cpp<double>(",
              self$expr(C[[1]]),
              "(",
              paste(vapply(as.list(C[-1]), self$expr, character(1)), collapse = ", "),
              "))"
            )
          }
        )
      } else {
        as.character(C)
      }
    },
    
    assign = function(var, value, type = "") {
      s <- paste0(var, " = ", value, ";")
      if (type != "") paste(type, s) else s
    },
    
    array = function(var, index) {
      paste0(var, "[", index, "]")
    },
    
    declare_array = function(var, length, type = "double") {
      paste0(type, " ", var, "[", length, "];")
    },
    
    block = function(statements, indent = "") {
      paste0(
        indent, "{", "\n",
        paste(paste0(indent, "  ", statements), collapse = "\n"),
        "\n", indent, "}"
      )
    },
    
    alias_statement = function(eq) {
      var <- eq[[2]]
      if (!is.name(var)) var <- as.name(var)
      self$assign(as.character(var), self$expr(eq[[3]]), "double")
    },
    
    bind_vars = function(S, name) {
      lapply(seq_along(S), function(i) {
        self$assign(S[[i]], self$array(paste0("__", name), i - 1), "double")
      })
    },
    
    substitutions = function() {
      lapply(self$alias, self$alias_statement)
    },
    
    switch_statement = function(var, values, indent) {
      l <- list()
      for (i in seq_along(values)) {
        n <- names(values)[i]
        if (is.null(n) || n == "") n <- i - 1
        l <- c(
          l,
          paste0("case ", n, ":"),
          paste0("  ", values[[i]]),
          "  break;"
        )
      }
      paste0("switch(", var, ") ", self$block(l, indent))
    },
    
    attached_functions = function() {
      fnames <- ls(attached.functions)
      if (length(fnames) == 0) return(NULL)
      
      lapply(fnames, function(n) {
        f <- attached.functions[[n]]
        if (is.null(f))
          stop("function ", n, " not attached.")
        paste0('cpp11::function ', n, '(attached_functions["', n, '"]);')
      })
    }
  )
)
