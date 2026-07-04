# ==============================================================================
# Gillespie: Gillespie SSA with optional compiled backend
# ==============================================================================
#
# This file implements a Simulator subclass that runs a Gillespie / SSA style
# stochastic simulation for a Compartmental model.
#
# `Gillespie` is the main class. It always builds the R step function exposed as
# `$model`; when `compile = TRUE`, it also builds a cpp11 backend in `$compiled`
# and uses that backend for simulation.
#
# `RGillespie` is retained as a compatibility wrapper for `compile = FALSE`.
# `CGillespie` is defined in CGillespie.R as a compatibility wrapper for
# `compile = TRUE`.
# ==============================================================================

#' R6 class implementing the Gillespie method
#'
#' This is a subclass of [Simulator], implementing the Gillespie stochastic
#' simulation algorithm for a [Compartmental] model.
#'
#' @docType class
#' @examples
#' # an SIR model
#' SIR <- Compartmental$new(S, I, R)
#' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
#' SIR$transition(I->R ~ gamma*I, name="recovery")
#' g <- Gillespie$new(SIR)
#' g$simulate(0:100, y0=c(S=1000, I=10, R=0), parms=c(beta=0.4, gamma=0.2))
#' @export
Gillespie <- R6Class(
  "Gillespie",
  inherit = Simulator,

  private = list(
    header = "#include \"cpp11.hpp\"
#include \"Rmath.h\"
#include <cmath>
",

    setup = "
cpp11::environment attached_functions(R_NilValue);

[[cpp11::register]]
void REpiSimSetup(cpp11::environment functions)
{
  attached_functions = functions;
}
",

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

    while (i < t.size() && time >= t[i]) {
      data(i, 0) = t[i];
      for (size_t j = 1; j <= y.size(); ++j) {
        data(i, j) = y[j-1];
      }
      ++i;
    }

    for (size_t j = 0; j < y.size(); ++j)
      y[j] = l[j+1];
  }

  return data;
}
",

    result = list(
      "cpp11::writable::doubles __result(__y.size() + 1);",
      "for (size_t i = 0; i < __y.size(); ++i)",
      "  __result[i+1] = __y[i];"
    ),

    middle = list(
      "if (__total == 0 || std::isnan(__total) || !std::isfinite(__total)) {",
      "  __result[0] = R_PosInf;",
      "  return __result;",
      "}",
      "double __p = unif_rand() * __total;",
      "size_t __transition = 0;",
      "for (; __rate[__transition] <= __p; ++__transition);"
    ),

    format.rate = function(T) {
      as.call(c(quote(c), lapply(unname(T), function(tr) tr$rate)))
    },

    format.substitution = function(S) {
      lapply(names(S), function(var) call("<-", as.name(var), S[[var]]))
    },

    format.var = function(S, name) {
      lapply(S, function(var) {
        call("<-", as.name(var), call("[[", as.name(name), match(var, S)))
      })
    },

    format.rate_cpp = function(T, cpp) {
      out <- list(cpp$declare_array("__rate", length(T)))
      for (i in seq_along(T)) {
        expr_i <- if (i > 1) {
          call("+", cpp$array("__rate", i - 2), T[[i]]$rate)
        } else {
          T[[i]]$rate
        }
        out <- c(out, cpp$assign(cpp$array("__rate", i - 1), cpp$expr(expr_i)))
      }
      c(out, cpp$assign("__total", cpp$array("__rate", length(T) - 1), "double"))
    },

    format.transition_cpp = function(transitions, cpp) {
      cases <- lapply(transitions, function(tr) {
        stmts <- list()

        if (!is.null(tr$from)) {
          idx_from <- which(private$compartments == tr$from)
          stmts <- c(stmts, cpp$assign(cpp$array("__result", idx_from), paste0(tr$from, " - 1")))
        }
        if (!is.null(tr$to)) {
          idx_to <- which(private$compartments == tr$to)
          stmts <- c(stmts, cpp$assign(cpp$array("__result", idx_to), paste0(tr$to, " + 1")))
        }

        stmts
      })

      cpp$switch_statement("__transition", unname(cases), indent = "  ")
    },

    build = function(model) {
      alias_names <- names(private$alias)
      extra_syms <- lapply(alias_names, as.name)
      names(extra_syms) <- alias_names

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

    build_step_cpp = function(model) {
      cpp <- CppConverter$new(private$compartments, private$parameters, private$alias)

      stmts <- c(
        cpp$bind_vars(private$compartments, "y"),
        cpp$bind_vars(private$parameters, "parms"),
        cpp$substitutions(),
        cpp$attached_functions(),
        private$result,
        private$format.rate_cpp(model$transitions, cpp),
        private$middle,
        cpp$assign("__result[0]", paste0(model$t, " + exp_rand()/__total")),
        private$format.transition_cpp(model$transitions, cpp),
        "return __result;"
      )

      paste0(
        "inline cpp11::doubles step(double ", model$t, ", cpp11::doubles __y, cpp11::doubles __parms) ",
        cpp$block(stmts, "")
      )
    },

    compile = function(model, r_model) {
      if (!requireNamespace("cpp11", quietly = TRUE)) {
        stop(
          "cpp11 is required to compile Gillespie. Use compile = FALSE or install cpp11",
          call. = FALSE
        )
      }

      program <- paste(
        private$header,
        private$setup,
        private$build_step_cpp(model),
        private$main,
        sep = "\n"
      )

      lib <- cpp11::cpp_source(code = program)
      simulate <- get("gillespie", inherits = TRUE)
      REpiSimSetup(attached.functions)

      list(
        program = program,
        lib = lib,
        simulate = simulate
      )
    },

    reorder_named_numeric = function(x, ord) {
      xx <- as.numeric(x[ord])
      names(xx) <- ord
      xx
    },

    check_initial_values = function(y) {
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
    },

    simulate_r = function(t, y, parms) {
      y <- y[private$compartments]
      parms <- parms[private$parameters]

      if (any(!is.finite(t)))
        stop("invalid time")
      if (length(t) < 2)
        stop("need at least two time points")
      if (is.null(names(y)) || any(names(y) == ""))
        stop("the initial values must be named")

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

        while (time > t[[i]]) {
          data[i, -1] <- y
          i <- i + 1
          if (i > length(t)) {
            return(data)
          }
        }

        y <- step[2:(1 + length(y))]
      }
    },

    simulate_cpp = function(t, y, parms) {
      y <- private$reorder_named_numeric(y, private$compartments)
      parms <- private$reorder_named_numeric(parms, private$parameters)

      data <- as.data.frame(private$.compiled$simulate(as.numeric(t), y, parms))
      colnames(data) <- c("time", names(y))
      data
    },

    .simulate = function(t, y, parms) {
      private$check_initial_values(y)

      if (!is.null(private$.compiled)) {
        private$simulate_cpp(t, y, parms)
      } else {
        private$simulate_r(t, y, parms)
      }
    },

    finalize = function() {
      lib <- private$.compiled$lib
      if (!is.null(lib) && !is.null(lib[["path"]])) {
        try(dyn.unload(lib[["path"]]), silent = TRUE)
      }
    }
  ),

  public = list(
    #' @description
    #' Construct a Gillespie simulator.
    #'
    #' @param model A [Compartmental] model.
    #' @param compile Logical. If `TRUE`, compile and use a cpp11 backend.
    initialize = function(model, compile = FALSE) {
      super$initialize(model, compile = compile)
    }
  )
)

#' R6 compatibility class for the R Gillespie backend
#'
#' `RGillespie$new(model)` is equivalent to `Gillespie$new(model, compile = FALSE)`.
#'
#' @docType class
#' @export
RGillespie <- R6Class(
  "RGillespie",
  inherit = Gillespie,

  public = list(
    initialize = function(model) {
      super$initialize(model, compile = FALSE)
    }
  )
)
