#' Declare a probability distribution family
#'
#' A distribution family is declared by its canonical density parameters and
#' the density function that evaluates them. Calling the family with parameter
#' values returns a `Distribution` object.
#'
#' If all canonical parameters are fixed, the returned object provides
#' `log.density(x)`. If one or more canonical parameters are missing, they are
#' computed from `mean` and any free parameters, and the returned object provides
#' `log.likelihood(x, mean, ...)`.
#'
#' @param canonical A named list of expressions defining each canonical density
#'   parameter in terms of `mean` and/or other canonical parameters.
#' @param density A density function accepting the canonical parameters and
#'   `log = TRUE`.
#' @return A distribution-family function.
#' @name Distribution
#' @export
Distribution <- function(canonical, density) {
  if (!is.list(canonical) || length(canonical) == 0)
    stop("canonical must be a non-empty named list")
  
  ns <- names(canonical)
  if (is.null(ns) || any(ns == ""))
    stop("canonical parameters must be named")
  if (!identical(ns, make.names(ns)))
    stop("canonical parameter names must be syntactic R names")
  
  if (!is.function(density))
    stop("density must be a function")
  
  canonical <- lapply(canonical, function(x) {
    if (is.null(x)) NULL else if (is.call(x) || is.name(x)) x else substitute(x)
  })
  
  family <- eval(parse(text = paste0("function(", paste(ns, collapse = ", "), ") NULL")))
  body(family) <- quote(distribution_instance(canonical, density, as.list(match.call())[-1]))
  
  class(family) <- c("DistributionFamily", class(family))
  attr(family, "canonical") <- canonical
  attr(family, "density") <- density
  family
}

#' @export
`$.DistributionFamily` <- function(x, name) {
  switch(
    name,
    canonical = attr(x, "canonical"),
    density = attr(x, "density"),
    new = x,
    stop("unknown distribution-family field: ", name)
  )
}

distribution_symbols <- function(expr) {
  if (is.null(expr)) return(character(0))
  setdiff(all.names(expr, functions = FALSE, unique = TRUE), "mean")
}

distribution_is_missing <- function(x) {
  missing(x) || is.null(x) || (length(x) == 1 && is.atomic(x) && is.na(x))
}

distribution_instance <- function(canonical, density, args) {
  cn <- names(canonical)
  supplied <- names(args)
  if (length(args) > length(cn))
    stop("too many distribution parameters")
  
  if (length(args) > 0) {
    if (is.null(supplied)) supplied <- rep("", length(args))
    unnamed <- which(supplied == "")
    if (length(unnamed) > 0)
      supplied[unnamed] <- cn[unnamed]
    names(args) <- supplied
  } else {
    supplied <- character(0)
  }
  
  bad <- setdiff(supplied, cn)
  if (length(bad) > 0)
    stop("unknown distribution parameter", if (length(bad) > 1) "s" else "",
         ": ", paste(bad, collapse = ", "))
  
  fixed <- list()
  missing_params <- character(0)
  
  for (p in cn) {
    if (!p %in% supplied || distribution_is_missing(args[[p]])) {
      missing_params <- c(missing_params, p)
    } else {
      fixed[[p]] <- args[[p]]
    }
  }
  
  if (length(missing_params) == 0) {
    log.density <- function(x) {
      do.call(density, c(list(x = x), fixed, list(log = TRUE)))
    }
    
    return(structure(
      list(
        canonical = canonical,
        density = density,
        values = fixed,
        parameters = character(0),
        par = character(0),
        log.density = log.density,
        log.likelihood = NULL
      ),
      class = "Distribution"
    ))
  }
  
  free <- character(0)
  derived <- list()
  available <- names(fixed)
  
  for (p in cn) {
    if (!p %in% missing_params) next
    
    expr <- canonical[[p]]
    symbols <- distribution_symbols(expr)
    
    if (is.null(expr) || (p != "mean" && identical(expr, as.name(p))) || p %in% symbols) {
      free <- union(free, p)
      available <- union(available, p)
      next
    }
    
    needed <- setdiff(symbols, available)
    if (length(needed) > 0) {
      unknown <- setdiff(needed, cn)
      if (length(unknown) > 0)
        stop("unknown symbol", if (length(unknown) > 1) "s" else "",
             " in canonical expression for ", p, ": ",
             paste(unknown, collapse = ", "))
      
      free <- union(free, needed)
      available <- union(available, needed)
    }
    
    derived[[p]] <- expr
    available <- union(available, p)
  }
  
  log.likelihood <- function(x, mean, ...) {
    extra <- list(...)
    extra_names <- names(extra)
    if (length(free) > 0 && (is.null(extra_names) || any(extra_names == "")))
      stop("distribution parameters must be named")
    
    missing_free <- setdiff(free, extra_names)
    if (length(missing_free) > 0)
      stop("missing distribution parameter", if (length(missing_free) > 1) "s" else "",
           ": ", paste(missing_free, collapse = ", "))
    
    unused <- setdiff(extra_names, free)
    if (length(unused) > 0)
      stop("unused distribution parameter", if (length(unused) > 1) "s" else "",
           ": ", paste(unused, collapse = ", "))
    
    env <- list2env(c(list(mean = mean), fixed, extra), parent = parent.frame())
    params <- vector("list", length(cn))
    names(params) <- cn
    
    for (p in cn) {
      if (!is.null(fixed[[p]])) {
        params[[p]] <- fixed[[p]]
      } else if (p %in% free) {
        params[[p]] <- get(p, envir = env, inherits = FALSE)
      } else {
        params[[p]] <- eval(derived[[p]], envir = env)
        assign(p, params[[p]], envir = env)
      }
    }
    
    sum(do.call(density, c(list(x = x), params, list(log = TRUE))))
  }
  
  structure(
    list(
      canonical = canonical,
      density = density,
      values = fixed,
      parameters = free,
      par = free,
      log.density = NULL,
      log.likelihood = log.likelihood
    ),
    class = "Distribution"
  )
}

#' The Poisson distribution
#'
#' @name Poisson
#' @export
Poisson <- Distribution(
  canonical = list(lambda = quote(mean)),
  density = dpois
)

#' The negative binomial distribution
#'
#' @name NBinom
#' @export
NBinom <- Distribution(
  canonical = list(
    size = quote(mean * prob / (1 - prob)),
    prob = quote(size / (mean + size))
  ),
  density = dnbinom
)

#' The normal distribution
#'
#' @name Normal
#' @export
Normal <- Distribution(
  canonical = list(
    mean = quote(mean),
    sd = quote(sd)
  ),
  density = dnorm
)

#' The uniform distribution
#'
#' @name Uniform
#' @export
Uniform <- Distribution(
  canonical = list(
    min = quote(2 * mean - max),
    max = quote(2 * mean + min)
  ),
  density = dunif
)

#' The exponential distribution
#'
#' @name Exponential
#' @export
Exponential <- Distribution(
  canonical = list(rate = quote(1 / mean)),
  density = dexp
)

#' The gamma distribution
#'
#' @name Gamma
#' @export
Gamma <- Distribution(
  canonical = list(
    shape = quote(mean / scale),
    scale = quote(mean / shape)
  ),
  density = dgamma
)

#' The beta distribution
#'
#' @name Beta
#' @export
Beta <- Distribution(
  canonical = list(
    shape1 = quote(mean / (1 - mean) * shape2),
    shape2 = quote((1 - mean) / mean * shape1)
  ),
  density = dbeta
)
