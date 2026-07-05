# ==============================================================================
# Stratification helpers
# ==============================================================================

strata_extract_index_args <- function(args, env) {
  ns <- names(args)
  if (is.null(ns)) ns <- rep("", length(args))
  idx <- which(ns == ".index")
  if (length(idx) > 1) stop(".index specified more than once")

  if (length(idx) == 0) {
    index_sets <- list()
  } else {
    index_sets <- strata_eval(args[[idx]], env)
    if (!is.list(index_sets) || is.null(names(index_sets)) || any(names(index_sets) == ""))
      stop(".index must be a named list")

    args <- args[-idx]
    ns <- ns[-idx]
  }

  inline_indices <- strata_formula_lhs_indices(args)
  inline_args <- which(ns %in% inline_indices)
  for (i in inline_args) {
    index_name <- ns[[i]]
    if (!is.null(index_sets[[index_name]]))
      stop("index ", index_name, " specified more than once")
    index_sets[[index_name]] <- strata_eval(args[[i]], env)
  }
  if (length(inline_args) > 0) args <- args[-inline_args]

  index_sets <- lapply(index_sets, as.character)
  list(args = args, index_sets = index_sets)
}

strata_formula_lhs_indices <- function(args) {
  indices <- character(0)
  for (arg in args) {
    if (!is.call(arg) || !identical(arg[[1]], as.name("~"))) next
    lhs <- arg[[2]]
    if (!strata_is_indexed_call(lhs)) next
    lhs_indices <- vapply(as.list(lhs)[-(1:2)], function(index) {
      if (!is.name(index)) stop("compartment indices must be symbols")
      as.character(index)
    }, character(1))
    indices <- union(indices, lhs_indices)
  }
  indices
}

strata_transition_endpoint_indices <- function(formula) {
  endpoints <- strata_transition_endpoints(formula)
  indices <- character(0)

  for (endpoint in endpoints) {
    if (!strata_is_indexed_call(endpoint)) next
    endpoint_indices <- vapply(as.list(endpoint)[-(1:2)], function(index) {
      if (!is.name(index)) stop("compartment indices must be symbols")
      as.character(index)
    }, character(1))
    indices <- union(indices, endpoint_indices)
  }

  indices
}

strata_extract_transition_index_args <- function(formula, args, env, global_index_sets = list()) {
  ns <- names(args)
  if (is.null(ns)) ns <- rep("", length(args))

  endpoint_indices <- strata_transition_endpoint_indices(formula)
  local_index_args <- which(ns %in% endpoint_indices)
  local_index_sets <- list()

  for (i in local_index_args) {
    index_name <- ns[[i]]
    if (!is.null(local_index_sets[[index_name]]))
      stop("index ", index_name, " specified more than once")
    local_index_sets[[index_name]] <- as.character(strata_eval(args[[i]], env))
  }

  if (length(local_index_args) > 0) args <- args[-local_index_args]

  index_sets <- global_index_sets
  for (index_name in names(local_index_sets)) {
    index_sets[[index_name]] <- local_index_sets[[index_name]]
  }

  for (index_name in endpoint_indices) {
    if (is.null(index_sets[[index_name]]))
      stop("index ", index_name, " must be specified for this transition")
  }

  list(args = args, index_sets = index_sets)
}

strata_is_indexed_call <- function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("[")) && length(expr) >= 3 &&
    is.name(expr[[2]])
}

strata_flat_name <- function(base, values) {
  make.names(paste(c(base, as.character(values)), collapse = "_"))
}

strata_flat_dimension_names <- function(name, dims) {
  grid <- expand.grid(dims, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  vapply(seq_len(nrow(grid)), function(row) {
    key <- as.character(unlist(grid[row, , drop = TRUE], use.names = FALSE))
    strata_flat_name(name, key)
  }, character(1))
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

strata_collect_parameter_dimensions <- function(formulas, index_sets, compartment_dimensions, env) {
  dimensions <- list()
  for (formula in formulas) {
    if (!is.call(formula) || !identical(formula[[1]], as.name("~")))
      stop("Invalid equation")

    lhs <- formula[[2]]
    bindings <- list()
    if (strata_is_indexed_call(lhs)) {
      base <- as.character(lhs[[2]])
      indices <- vapply(as.list(lhs)[-(1:2)], function(index) {
        if (!is.name(index)) stop("compartment indices must be symbols")
        as.character(index)
      }, character(1))
      dims <- compartment_dimensions[[base]]
      for (i in seq_along(indices))
        bindings[[indices[[i]]]] <- dims[[i]]
    }

    dimensions <- strata_collect_parameter_dimensions_expr(
      formula[[3]], bindings, index_sets, compartment_dimensions, env, dimensions
    )
  }
  dimensions
}

strata_collect_transition_compartment_dimensions <- function(formula, index_sets, dimensions = list()) {
  endpoints <- strata_transition_endpoints(formula)

  for (endpoint in endpoints) {
    if (!strata_is_indexed_call(endpoint)) next

    base <- as.character(endpoint[[2]])
    indices <- as.list(endpoint)[-(1:2)]
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

strata_collect_transition_parameter_dimensions <- function(formula, index_sets,
                                                           compartment_dimensions,
                                                           env,
                                                           dimensions = list()) {
  pieces <- strata_parse_transition_formula(formula)
  bindings <- list()

  for (endpoint in list(pieces$to, pieces$from)) {
    if (!strata_is_indexed_call(endpoint)) next

    base <- as.character(endpoint[[2]])
    indices <- vapply(as.list(endpoint)[-(1:2)], function(index) {
      if (!is.name(index)) stop("compartment indices must be symbols")
      as.character(index)
    }, character(1))
    dims <- compartment_dimensions[[base]]
    for (i in seq_along(indices))
      bindings[[indices[[i]]]] <- dims[[i]]
  }

  strata_collect_parameter_dimensions_expr(
    pieces$rate, bindings, index_sets, compartment_dimensions, env, dimensions
  )
}

strata_transition_endpoints <- function(formula) {
  pieces <- strata_parse_transition_formula(formula)
  Filter(Negate(is.null), list(pieces$to, pieces$from))
}

strata_parse_transition_formula <- function(formula) {
  if (!is.call(formula) || !identical(formula[[1]], as.name("<-")))
    stop("invalid transition")

  to <- strata_parse_transition_side(formula[[2]])
  from <- strata_parse_transition_side(formula[[3]])

  rate <- NULL
  if (!is.null(to$rate)) rate <- to$rate
  if (!is.null(from$rate)) {
    if (!is.null(rate)) stop("invalid transition: rate specified twice")
    rate <- from$rate
  }

  list(to = to$compartment, from = from$compartment, rate = rate)
}

strata_parse_transition_side <- function(side) {
  if (is.null(side)) return(list(compartment = NULL, rate = NULL))
  if (is.name(side) && identical(as.character(side), "NULL"))
    return(list(compartment = NULL, rate = NULL))
  if (is.call(side) && identical(side[[1]], as.name("~"))) {
    return(list(compartment = side[[2]], rate = side[[3]]))
  }
  list(compartment = side, rate = NULL)
}

strata_collect_parameter_dimensions_expr <- function(expr, bindings, index_sets,
                                                     compartment_dimensions, env,
                                                     dimensions) {
  if (is.name(expr)) {
    return(dimensions)
  }

  if (!is.call(expr)) return(dimensions)

  fn <- as.character(expr[[1]])
  if (fn %in% c("Sum", "sum", "Prod", "prod") && strata_has_index_bindings(expr)) {
    args <- as.list(expr)[-1]
    index_args <- args[-1]
    index_names <- names(index_args)
    local_bindings <- bindings
    for (i in seq_along(index_args)) {
      index_name <- index_names[[i]]
      values <- index_sets[[index_name]]
      if (is.null(values))
        values <- strata_eval(index_args[[i]], env)
      local_bindings[[index_name]] <- as.character(values)
    }
    return(strata_collect_parameter_dimensions_expr(
      args[[1]], local_bindings, index_sets, compartment_dimensions, env, dimensions
    ))
  }

  if (strata_is_indexed_call(expr)) {
    base <- as.character(expr[[2]])
    if (is.null(compartment_dimensions[[base]])) {
      indices <- as.list(expr)[-(1:2)]
      dims <- lapply(indices, function(index) {
        if (!is.name(index)) return(NULL)
        index_name <- as.character(index)
        values <- bindings[[index_name]]
        if (is.null(values)) values <- index_sets[[index_name]]
        if (is.null(values)) return(NULL)
        as.character(values)
      })

      if (all(!vapply(dims, is.null, logical(1)))) {
        old <- dimensions[[base]]
        if (!is.null(old) && !identical(old, dims))
          stop("conflicting dimensions for parameter ", base)
        dimensions[[base]] <- dims
      }
    }
    return(dimensions)
  }

  parts <- as.list(expr)
  for (i in seq_along(parts)[-1]) {
    dimensions <- strata_collect_parameter_dimensions_expr(
      parts[[i]], bindings, index_sets, compartment_dimensions, env, dimensions
    )
  }
  dimensions
}

strata_expand_model_formula <- function(formula, index_sets, compartment_dimensions, env,
                                        parameter_dimensions = list(),
                                        index_mode = c("name", "position")) {
  index_mode <- match.arg(index_mode)
  if (!is.call(formula) || !identical(formula[[1]], as.name("~")))
    stop("Invalid equation")

  lhs <- formula[[2]]
  rhs <- formula[[3]]

  if (!strata_is_indexed_call(lhs)) {
    return(list(call("~", lhs, strata_expand_expr(
      rhs, list(), index_sets, compartment_dimensions, env, parameter_dimensions, index_mode
    ))))
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
    positions <- vapply(seq_along(dims), function(i) match(values[[i]], dims[[i]]), integer(1))
    names(positions) <- indices
    current_position <- if (length(positions) == 1) positions[[1]] else NULL
    call(
      "~",
      as.name(strata_flat_name(base, values)),
      strata_expand_expr(
        rhs,
        strata_make_bindings(indices, values, positions),
        index_sets,
        compartment_dimensions,
        env,
        parameter_dimensions,
        index_mode,
        current_position
      )
    )
  })
}

strata_expand_transition_formula <- function(formula, index_sets, compartment_dimensions, env,
                                             parameter_dimensions = list(),
                                             index_mode = c("name", "position")) {
  index_mode <- match.arg(index_mode)

  if (length(index_sets) == 0)
    return(list(list(formula = formula, suffix = NULL)))

  pieces <- strata_parse_transition_formula(formula)
  endpoint_indices <- list()
  for (endpoint in list(pieces$to, pieces$from)) {
    if (!strata_is_indexed_call(endpoint)) next
    indices <- vapply(as.list(endpoint)[-(1:2)], function(index) {
      if (!is.name(index)) stop("compartment indices must be symbols")
      as.character(index)
    }, character(1))
    for (index in indices) endpoint_indices[[index]] <- index_sets[[index]]
  }

  if (length(endpoint_indices) == 0)
    return(list(list(
      formula = strata_expand_expr(
        formula, list(), index_sets, compartment_dimensions, env,
        parameter_dimensions, index_mode
      ),
      suffix = NULL
    )))

  index_names <- names(endpoint_indices)
  dims <- unname(endpoint_indices)
  grid <- expand.grid(dims, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  lapply(seq_len(nrow(grid)), function(row) {
    values <- as.character(unlist(grid[row, , drop = TRUE], use.names = FALSE))
    positions <- vapply(seq_along(dims), function(i) match(values[[i]], dims[[i]]), integer(1))
    names(values) <- index_names
    suffix <- paste(values, collapse = "_")
    current_position <- if (length(positions) == 1) positions[[1]] else NULL
    bindings <- strata_make_bindings(index_names, unname(values), positions)
    expanded_to <- if (is.null(pieces$to)) NULL else strata_expand_expr(
      pieces$to, bindings, index_sets, compartment_dimensions, env,
      parameter_dimensions, index_mode, current_position
    )
    expanded_from <- if (is.null(pieces$from)) NULL else strata_expand_expr(
      pieces$from, bindings, index_sets, compartment_dimensions, env,
      parameter_dimensions, index_mode, current_position
    )
    expanded_rate <- if (is.null(pieces$rate)) NULL else strata_expand_expr(
      pieces$rate, bindings, index_sets, compartment_dimensions, env,
      parameter_dimensions, index_mode, current_position
    )

    list(
      formula = call("<-", expanded_to, call("~", expanded_from, expanded_rate)),
      suffix = suffix
    )
  })
}

strata_expand_expr <- function(expr, bindings, index_sets, compartment_dimensions, env,
                               parameter_dimensions = list(),
                               index_mode = c("name", "position"),
                               current_position = NULL) {
  index_mode <- match.arg(index_mode)
  if (is.name(expr)) {
    name <- as.character(expr)
    if (!is.null(bindings[[name]]))
      stop("bare index symbol '", name, "' is not allowed outside indexed expressions")
    if (index_mode == "position" && !is.null(current_position) &&
        !is.null(parameter_dimensions[[name]]) && length(parameter_dimensions[[name]]) == 1) {
      values <- parameter_dimensions[[name]][[1]]
      return(as.name(strata_flat_name(name, values[[current_position]])))
    }
    return(expr)
  }

  if (!is.call(expr)) return(expr)

  fn <- as.character(expr[[1]])
  if (fn %in% c("Sum", "sum") && strata_has_index_bindings(expr))
    return(strata_expand_reduction(
      expr, "+", bindings, index_sets, compartment_dimensions, env,
      parameter_dimensions, index_mode, current_position
    ))
  if (fn %in% c("Prod", "prod") && strata_has_index_bindings(expr))
    return(strata_expand_reduction(
      expr, "*", bindings, index_sets, compartment_dimensions, env,
      parameter_dimensions, index_mode, current_position
    ))

  if (strata_is_indexed_call(expr)) {
    base <- as.character(expr[[2]])
    indices <- as.list(expr)[-(1:2)]
    values <- lapply(indices, function(index) {
      if (is.name(index)) {
        index_name <- as.character(index)
        binding <- bindings[[index_name]]
        if (is.null(binding)) stop("unbound index: ", index_name)
        if (index_mode == "position") binding$position else binding$value
      } else {
        strata_expand_expr(index, bindings, index_sets, compartment_dimensions, env, parameter_dimensions, index_mode, current_position)
      }
    })

    if (!is.null(compartment_dimensions[[base]])) {
      if (length(values) != length(compartment_dimensions[[base]]))
        stop("wrong number of indices for compartment ", base)
      flat_values <- lapply(indices, function(index) {
        if (!is.name(index)) stop("compartment indices must be symbols")
        binding <- bindings[[as.character(index)]]
        if (is.null(binding)) stop("unbound index: ", as.character(index))
        binding$value
      })
      return(as.name(strata_flat_name(base, unlist(flat_values, use.names = FALSE))))
    }

    if (!is.null(parameter_dimensions[[base]])) {
      if (length(values) != length(parameter_dimensions[[base]]))
        stop("wrong number of indices for parameter ", base)
      flat_values <- lapply(indices, function(index) {
        if (is.name(index)) {
          binding <- bindings[[as.character(index)]]
          if (is.null(binding)) stop("unbound index: ", as.character(index))
          return(binding$value)
        }

        value <- strata_expand_expr(index, bindings, index_sets, compartment_dimensions, env, parameter_dimensions, "name", current_position)
        if (!is.character(value) && !is.numeric(value))
          stop("parameter indices must be symbols or scalar values")
        as.character(value)
      })
      return(as.name(strata_flat_name(base, unlist(flat_values, use.names = FALSE))))
    }

    return(as.call(c(list(as.name("["), as.name(base)), unname(values))))
  }

  parts <- as.list(expr)
  for (i in seq_along(parts)[-1])
    parts[[i]] <- strata_expand_expr(parts[[i]], bindings, index_sets, compartment_dimensions, env, parameter_dimensions, index_mode, current_position)
  as.call(parts)
}

strata_has_index_bindings <- function(expr) {
  arg_names <- names(as.list(expr))
  !is.null(arg_names) && any(arg_names[-1] != "")
}

strata_expand_reduction <- function(expr, op, bindings, index_sets, compartment_dimensions, env,
                                    parameter_dimensions = list(),
                                    index_mode = c("name", "position"),
                                    current_position = NULL) {
  index_mode <- match.arg(index_mode)
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
    positions <- vapply(seq_along(dims), function(i) match(values[[i]], dims[[i]]), integer(1))
    strata_expand_expr(
      body,
      c(bindings, strata_make_bindings(index_names, values, positions)),
      index_sets,
      compartment_dimensions,
      env,
      parameter_dimensions,
      index_mode,
      current_position
    )
  })

  if (length(terms) == 0)
    return(if (op == "+") 0 else 1)
  Reduce(function(a, b) call(op, a, b), terms)
}

strata_position_or_scalar_call <- function(expr, position) {
  call(
    "[[",
    expr,
    call("min", call("length", expr), as.integer(position))
  )
}

strata_make_bindings <- function(index_names, values, positions) {
  bindings <- lapply(seq_along(index_names), function(i) {
    list(value = values[[i]], position = positions[[i]])
  })
  names(bindings) <- index_names
  bindings
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
