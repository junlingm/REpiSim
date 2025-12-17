# ==============================================================================
# CGillespie: Gillespie SSA (cpp11 backend)
# ==============================================================================
#
# This file implements a Simulator subclass that runs a Gillespie / SSA style
# stochastic simulation for a Compartmental model using C++ via cpp11.
#
# High-level approach:
#   1) Compile model expressions (rates + substitutions) into an inline C++
#      function step(time, y, parms) that advances the system by one event.
#   2) Expose a registered C++ function gillespie(t, y0, parms) that repeatedly
#      calls step() and records the state on a requested time grid `t`.
#   3) Compile and load the generated C++ code at runtime using cpp11::cpp_source().
#
# NOTE (security):
#   This implementation calls user-attached R functions from C++ by holding them
#   as cpp11::function objects. In hosted environments (e.g., Shiny), consider
#   using a restricted/locked function registry and expression validation.
# ==============================================================================

#' R6 class implementing the Gillespie method using cpp11
#'
#' This is a subclass of [Simulator], implementing the Gillespie method to
#' simulate a [Compartmental] model using C++ via **cpp11**. This class is only
#' available if the `cpp11` package is installed.
#'
#' @docType class
#' @examples
#' # an SIR model
#' SIR <- Compartmental$new(S, I, R)
#' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
#' SIR$transition(I->R ~ gamma*I, name="recovery")
#' g <- CGillespie$new(SIR)
#' g$simulate(0:100, y0=c(S=1000, I=10, R=0), parms=c(beta=0.4, gamma=0.2))
#' @export
CGillespie <- R6Class(
  "CGillespie",
  inherit = Simulator,
  
  private = list(
    # --------------------------------------------------------------------------
    # Runtime-compiled artifacts
    # --------------------------------------------------------------------------
    
    # Full generated C++ program text (for debugging / inspection)
    .program = "",
    
    # Dynamic library returned by cpp11::cpp_source()
    lib = NULL,
    
    # Handle to compiled gillespie() function (an R function calling into C++)
    gillespie = NULL,
    
    # --------------------------------------------------------------------------
    # C++ template code
    # --------------------------------------------------------------------------
    
    header = "#include \"cpp11.hpp\"
#include \"Rmath.h\"
#include <cmath>
",
    
    # Bridge to pass the attached.functions environment into C++
    setup = "
cpp11::environment attached_functions(R_NilValue);

[[cpp11::register]]
void REpiSimSetup(cpp11::environment functions)
{
  attached_functions = functions;
}
",
    
    # Main loop: simulate on a time grid t, calling step() repeatedly
    main = "[[cpp11::register]]
cpp11::doubles_matrix<> gillespie(cpp11::doubles t, cpp11::doubles y0, cpp11::doubles parms) {
  cpp11::writable::doubles y(y0);
  cpp11::writable::doubles_matrix<> data(t.size(), y.size() + 1);

  if (!y.named())
    cpp11::stop(\"the initial values must be named\");

  double time = t[0];
  size_t i = 0;

  while (i < t.size()) {
    auto l = step(time, y, parms);
    time = l[0];

    // record state at all requested grid points that have been reached
    while (i < t.size() && time >= t[i]) {
      data(i, 0) = t[i];
      for (size_t j = 1; j <= y.size(); ++j) {
        data(i, j) = y[j-1];
      }
      ++i;
    }

    // update state
    for (size_t j = 0; j < y.size(); ++j)
      y[j] = l[j+1];
  }

  return data;
}
",

# --------------------------------------------------------------------------
# Code fragments used to build step()
# --------------------------------------------------------------------------

# Allocate result vector and copy y into it (shifted by 1 for time)
result = list(
  "cpp11::writable::doubles __result(__y.size() + 1);",
  "for (size_t i = 0; i < __y.size(); ++i)",
  "  __result[i+1] = __y[i];"
),

# Sample next transition, given cumulative rates in __rate[ ]
middle = list(
  "if (__total == 0 || std::isnan(__total) || !std::isfinite(__total)) {",
  "  __result[0] = R_PosInf;",
  "  return __result;",
  "}",
  "double __p = unif_rand() * __total;",
  "size_t __transition = 0;",
  "for (; __rate[__transition] <= __p; ++__transition);"
),

# --------------------------------------------------------------------------
# Small utilities for code generation
# --------------------------------------------------------------------------

# Remove bracket wrapper from a bracketed expression
remove.bracket = function(x) {
  if (is.call(x) && identical(x[[1]], as.name("(")))
    private$remove.bracket(x[[2]])
  else x
},

# Convert a simple R expression AST into a C++ expression string.
#
# Supported operators: + - * / ^ (), indexing, and function calls.
# Function calls are emitted as cpp11::as_cpp<double>(fn(args...)) so that
# both C++ functions and cpp11::function wrappers can be used uniformly.
cpp = function(C) {
  if (is.call(C)) {
    # normalize R's [[ to [ for our purposes
    if (identical(C[[1]], as.name("[["))) C[[1]] <- as.name("[")
    
    op <- as.character(C[[1]])
    
    switch(
      op,
      `+` = if (length(C) == 2) private$cpp(C[[2]]) else
        paste0(private$cpp(C[[2]]), " + ", private$cpp(C[[3]])),
      `-` = if (length(C) == 2) {
        paste0("-", private$cpp(C[[2]]))
      } else {
        paste0(private$cpp(C[[2]]), " - ", private$cpp(C[[3]]))
      },
      `*` = paste0(private$cpp(C[[2]]), " * ", private$cpp(C[[3]])),
      `/` = paste0(private$cpp(C[[2]]), " / ", private$cpp(C[[3]])),
      `^` = paste0(
        "pow(",
        private$cpp(private$remove.bracket(C[[2]])),
        ", ",
        private$cpp(private$remove.bracket(C[[3]])),
        ")"
      ),
      `(` = paste0("(", private$cpp(C[[2]]), ")"),
      `[` = if (length(C) > 3) {
        # support nested indexing like x[i][j] by recursively emitting calls
        private$cpp(C[-1])
      } else {
        paste0(private$cpp(C[[2]]), "[", private$cpp(C[[3]]), "]")
      },
      {
        # generic function call
        paste0(
          "cpp11::as_cpp<double>(",
          private$cpp(C[[1]]),
          "(",
          paste(vapply(as.list(C[-1]), private$cpp, character(1)), collapse = ", "),
          "))"
        )
      }
    )
  } else {
    as.character(C)
  }
},

# Generate a C++ assignment statement; optionally declare type.
assign = function(var, value, type = "") {
  s <- paste0(var, " = ", value, ";")
  if (type != "") paste(type, s) else s
},

# Generate C++ array access (by index)
array = function(var, index) {
  paste0(var, "[", index, "]")
},

# Declare a fixed-size C++ array
declare_array = function(var, length, type = "double") {
  paste0(type, " ", var, "[", length, "];")
},

# Format transition rate calculation into cumulative rates __rate[ ] and __total.
#
# We store cumulative rates to make transition selection O(M) with one scan.
format.rate = function(T) {
  out <- list(private$declare_array("__rate", length(T)))
  for (i in seq_along(T)) {
    expr_i <- if (i > 1) {
      call("+", private$array("__rate", i - 2), T[[i]]$rate)
    } else {
      T[[i]]$rate
    }
    out <- c(out, private$assign(private$array("__rate", i - 1), private$cpp(expr_i)))
  }
  c(out, private$assign("__total", private$array("__rate", length(T) - 1), "double"))
},

# Build a C++ block with indentation.
block = function(statements, indent = "") {
  paste0(
    indent, "{", "\n",
    paste(paste0(indent, "  ", statements), collapse = "\n"),
    "\n", indent, "}"
  )
},

# Format a model equation (alias assignment) as a C++ statement.
#
# eq is an R call shaped like:  (<-) name expr
format.alias_stmt = function(eq) {
  var <- eq[[2]]
  if (!is.name(var)) var <- as.name(var)
  private$assign(as.character(var), private$cpp(eq[[3]]), "double")
},

# Bind y and parms into C++ locals:
#   double S = __y[0]; ...
#   double beta = __parms[0]; ...
format.var_cpp = function(S, name) {
  lapply(seq_along(S), function(i) {
    private$assign(S[[i]], private$array(paste0("__", name), i - 1), "double")
  })
},

# Format substitutions (aliases) into C++ locals.
#
# We use the already dependency-sorted alias list from Simulator (private$alias),
# and generate: double N = S + I + R; etc.
format.substitution_cpp = function() {
  lapply(private$alias, private$format.alias_stmt)
},

# Format switch statement
format.switch = function(var, values, indent) {
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
  paste0("switch(", var, ") ", private$block(l, indent))
},

# Format state updates for each transition as cases in a switch statement.
#
# __result already contains the current y. We update the relevant entries:
#   from: __result[idx] = from - 1;
#   to  : __result[idx] = to + 1;
format.transition = function(transitions) {
  cases <- lapply(transitions, function(tr) {
    stmts <- list()
    
    if (!is.null(tr$from)) {
      idx_from <- which(private$compartments == tr$from)
      stmts <- c(stmts, private$assign(private$array("__result", idx_from), paste0(tr$from, " - 1")))
    }
    if (!is.null(tr$to)) {
      idx_to <- which(private$compartments == tr$to)
      stmts <- c(stmts, private$assign(private$array("__result", idx_to), paste0(tr$to, " + 1")))
    }
    
    stmts
  })
  
  private$format.switch("__transition", unname(cases), indent = "  ")
},

# Declare attached functions (custom R functions callable from C++).
#
# Each attached function becomes: cpp11::function fname(attached_functions["fname"]);
attach.functions = function() {
  fnames <- ls(attached.functions)
  if (length(fnames) == 0) return(NULL)
  
  lapply(fnames, function(n) {
    f <- attached.functions[[n]]
    if (is.null(f))
      stop("function ", n, " not attached.")
    paste0('cpp11::function ', n, '(attached_functions["', n, '"]);')
  })
},

# Build the C++ inline step() function.
build = function(model) {
  stmts <- c(
    # bind states and parameters
    private$format.var_cpp(private$compartments, "y"),
    private$format.var_cpp(private$parameters, "parms"),
    
    # compute substitutions (aliases)
    private$format.substitution_cpp(),
    
    # bring attached R functions into scope (if any)
    private$attach.functions(),
    
    # allocate result
    private$result,
    
    # compute cumulative rates + total
    private$format.rate(model$transitions),
    
    # sample event
    private$middle,
    
    # sample next time increment (Exp(total))
    private$assign("__result[0]", paste0(model$t, " + exp_rand()/__total")),
    
    # apply transition (±1) and return result
    private$format.transition(model$transitions),
    "return __result;"
  )
  
  paste0(
    "inline cpp11::doubles step(double ", model$t, ", cpp11::doubles __y, cpp11::doubles __parms) ",
    private$block(stmts, "")
  )
},

# Reorder numeric vector/list by desired name order (stable helper).
reorder_named_numeric = function(x, ord) {
  xx <- as.numeric(x[ord])
  names(xx) <- ord
  xx
},

# Perform the simulation by calling the compiled C++ gillespie().
.simulate = function(t, y, parms) {
  # validate integer, nonnegative initial conditions
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
  
  y <- private$reorder_named_numeric(y, private$compartments)
  parms <- private$reorder_named_numeric(parms, private$parameters)
  
  data <- as.data.frame(private$gillespie(as.numeric(t), y, parms))
  colnames(data) <- c("time", names(y))
  data
},

# Destructor: unload the compiled dynamic library (best effort).
finalize = function() {
  if (!is.null(private$lib) && !is.null(private$lib[["path"]])) {
    try(dyn.unload(private$lib[["path"]]), silent = TRUE)
  }
}
  ),

public = list(
  #' @description
  #' Construct a C++ Gillespie simulator.
  #'
  #' @param model A [Compartmental] model.
  #' @details
  #' `CGillespie` is available only if the **cpp11** package is installed.
  initialize = function(model) {
    if (!requireNamespace("cpp11", quietly = TRUE)) {
      stop(
        "cpp11 is required to use CGillespie. Either install the package, or use RGillespie instead",
        call. = FALSE
      )
    }
    
    super$initialize(model)
    
    # Compose full C++ program.
    private$.program <- paste(
      private$header,
      private$setup,
      self$model,   # compiled step()
      private$main, # compiled gillespie()
      sep = "\n"
    )
    
    # Compile and load.
    private$lib <- cpp11::cpp_source(code = private$.program)
    
    # cpp_source registers gillespie() into the calling environment.
    private$gillespie <- get("gillespie", inherits = TRUE)
    
    # Transfer attached.functions env to C++.
    REpiSimSetup(attached.functions)
  }
)
)
