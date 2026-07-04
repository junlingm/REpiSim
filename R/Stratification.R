# ==============================================================================
# Stratification helpers
# ==============================================================================

strata_extract_index_args <- function(args, env) {
  ns <- names(args)
  if (is.null(ns)) ns <- rep("", length(args))
  idx <- which(ns == ".index")
  if (length(idx) == 0) return(list(args = args, index_sets = list()))
  if (length(idx) > 1) stop(".index specified more than once")

  index_sets <- strata_eval(args[[idx]], env)
  if (!is.list(index_sets) || is.null(names(index_sets)) || any(names(index_sets) == ""))
    stop(".index must be a named list")

  index_sets <- lapply(index_sets, as.character)
  args <- args[-idx]
  list(args = args, index_sets = index_sets)
}

strata_is_indexed_call <- function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("[")) && length(expr) >= 3 &&
    is.name(expr[[2]])
}

strata_flat_name <- function(base, values) {
  make.names(paste(c(base, as.character(values)), collapse = "_"))
}

strata_collect_compartment_dimensions <- function(formulas, index_sets) {
  dimensions <- list()

  for (formula in formulas) {
    if (!is.call(formula) || !identical(formula[[1]], as.name("~")))
      stop("Invalid equation")

    lhs <- formula[[2]]
    if (!strata_is_indexed_call(lhs)) next

    base <- as.character(lhs[[2]])
    indices <- as.list(lhs)[-(1:2)]
    dims <- lapply(indices, function(index) {
      if (!is.name(index)) stop("compartment indices must be symbols")
      index_name <- as.character(index)
      values <- index_sets[[index_name]]
      if (is.null(values)) stop("unknown index: ", index_name)
      values
    })

    old <- dimensions[[base]]
    if (!is.null(old) && !identical(old, dims))
      stop("conflicting dimensions for compartment ", base)
    dimensions[[base]] <- dims
  }

  dimensions
}

strata_expand_model_formula <- function(formula, index_sets, compartment_dimensions, env) {
  if (!is.call(formula) || !identical(formula[[1]], as.name("~")))
    stop("Invalid equation")

  lhs <- formula[[2]]
  rhs <- formula[[3]]

  if (!strata_is_indexed_call(lhs)) {
    return(list(call("~", lhs, strata_expand_expr(rhs, list(), index_sets, compartment_dimensions, env))))
  }

  base <- as.character(lhs[[2]])
  indices <- vapply(as.list(lhs)[-(1:2)], function(index) {
    if (!is.name(index)) stop("compartment indices must be symbols")
    as.character(index)
  }, character(1))

  dims <- compartment_dimensions[[base]]
  if (is.null(dims)) stop("unknown indexed compartment ", base)
  if (length(indices) != length(dims))
    stop("wrong number of indices for compartment ", base)

  grid <- expand.grid(dims, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  lapply(seq_len(nrow(grid)), function(row) {
    values <- as.character(unlist(grid[row, , drop = TRUE], use.names = FALSE))
    names(values) <- indices
    call(
      "~",
      as.name(strata_flat_name(base, values)),
      strata_expand_expr(rhs, as.list(values), index_sets, compartment_dimensions, env)
    )
  })
}

strata_expand_expr <- function(expr, bindings, index_sets, compartment_dimensions, env) {
  if (is.name(expr)) {
    name <- as.character(expr)
    if (!is.null(bindings[[name]]))
      stop("bare index symbol '", name, "' is not allowed outside indexed expressions")
    return(expr)
  }

  if (!is.call(expr)) return(expr)

  fn <- as.character(expr[[1]])
  if (fn %in% c("Sum", "sum") && strata_has_index_bindings(expr))
    return(strata_expand_reduction(expr, "+", bindings, index_sets, compartment_dimensions, env))
  if (fn %in% c("Prod", "prod") && strata_has_index_bindings(expr))
    return(strata_expand_reduction(expr, "*", bindings, index_sets, compartment_dimensions, env))

  if (strata_is_indexed_call(expr)) {
    base <- as.character(expr[[2]])
    indices <- as.list(expr)[-(1:2)]
    values <- lapply(indices, function(index) {
      if (is.name(index)) {
        index_name <- as.character(index)
        value <- bindings[[index_name]]
        if (is.null(value)) stop("unbound index: ", index_name)
        value
      } else {
        strata_expand_expr(index, bindings, index_sets, compartment_dimensions, env)
      }
    })

    if (!is.null(compartment_dimensions[[base]])) {
      if (length(values) != length(compartment_dimensions[[base]]))
        stop("wrong number of indices for compartment ", base)
      return(as.name(strata_flat_name(base, unlist(values, use.names = FALSE))))
    }

    return(as.call(c(list(as.name("["), as.name(base)), unname(values))))
  }

  parts <- as.list(expr)
  for (i in seq_along(parts)[-1])
    parts[[i]] <- strata_expand_expr(parts[[i]], bindings, index_sets, compartment_dimensions, env)
  as.call(parts)
}

strata_has_index_bindings <- function(expr) {
  arg_names <- names(as.list(expr))
  !is.null(arg_names) && any(arg_names[-1] != "")
}

strata_expand_reduction <- function(expr, op, bindings, index_sets, compartment_dimensions, env) {
  args <- as.list(expr)[-1]
  ns <- names(args)
  if (length(args) == 0 || is.null(ns) || ns[[1]] != "")
    stop(expr[[1]], "() requires an expression followed by named index bindings")

  body <- args[[1]]
  index_args <- args[-1]
  index_names <- names(index_args)
  if (is.null(index_names) || any(index_names == ""))
    stop(expr[[1]], "() index bindings must be named")

  dims <- lapply(seq_along(index_args), function(i) {
    index_name <- index_names[[i]]
    values <- index_sets[[index_name]]
    if (is.null(values))
      values <- strata_eval(index_args[[i]], env)
    as.character(values)
  })

  grid <- expand.grid(dims, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  terms <- lapply(seq_len(nrow(grid)), function(row) {
    values <- as.character(unlist(grid[row, , drop = TRUE], use.names = FALSE))
    names(values) <- index_names
    strata_expand_expr(body, c(bindings, as.list(values)), index_sets, compartment_dimensions, env)
  })

  if (length(terms) == 0)
    return(if (op == "+") 0 else 1)
  Reduce(function(a, b) call(op, a, b), terms)
}

strata_eval <- function(expr, env) {
  value <- try(eval(expr, envir = env), silent = TRUE)
  if (!inherits(value, "try-error")) return(value)

  frames <- rev(sys.frames())
  for (frame in frames) {
    value <- try(eval(expr, envir = frame), silent = TRUE)
    if (!inherits(value, "try-error")) return(value)
  }

  eval(expr, envir = env)
}
