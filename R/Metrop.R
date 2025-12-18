# ==============================================================================
# MCMC calibration via mcmc::metrop
# ==============================================================================
#
# This file defines the Metrop calibrator, which runs a Metropolis sampler using
# mcmc::metrop and summarizes posterior samples using coda.
#
# PATCH NOTE (2025-12-17):
#   - Compatible with Calibrator storing private$.data as a data.frame even when
#     only a single observation column is mapped.
#   - Handles the case where the simulator returns a vector for a single series
#     (e.g., cumulative data converted to incidence via diff()).
# ==============================================================================

#' MCMC calibrator using `mcmc::metrop`
#'
#' @details
#' `Metrop` is a Bayesian calibrator that samples from the posterior distribution
#' using the Metropolis algorithm implemented in [mcmc::metrop()].
#'
#' The log-posterior is computed as:
#' \deqn{\log p(\theta \mid y) = \log p(y \mid \theta) + \log p(\theta)}
#' where:
#' - \eqn{p(y \mid \theta)} is specified by a `Distribution` likelihood object,
#' - \eqn{p(\theta)} is specified by a named list of `Distribution` prior objects.
#'
#' @name Metrop
#' @docType class
#' @export
Metrop <- R6::R6Class(
  "Metrop",
  inherit = Bayesian,
  
  private = list(
    #' @field .likelihood Distribution object defining the likelihood
    .likelihood = NULL,
    
    #' Compute log-posterior for a proposed parameter vector.
    #'
    #' @param pars Numeric vector proposed by `mcmc::metrop`.
    #' @param formula Named list of Expression objects to compute derived quantities.
    #' @param fixed Named list of fixed values.
    #' @param priors Named list of prior Distribution objects (names match pars).
    #' @param ... Forwarded to the underlying simulator.
    objective = function(pars, formula, fixed, priors, ...) {
      names(pars) <- names(priors)
      all <- c(as.list(pars), fixed)
      
      # Log prior
      lp <- sum(vapply(names(priors), function(n) priors[[n]]$log.density(all[[n]]), numeric(1)))
      if (!is.finite(lp)) return(-Inf)
      
      # Derived quantities (dependency order ensured upstream)
      for (n in names(formula)) {
        f <- formula[[n]]
        all[[n]] <- eval(f$expr, envir = all)
      }
      
      # Likelihood parameters (if used by the likelihood Distribution)
      pars_l <- all[private$.par.likelihood]
      
      # Simulated model outputs aligned with mapping
      x <- private$simulate(all, NULL, NULL, ...)
      
      # Observed data (always a data.frame in patched Calibrator)
      yobs <- private$.data
      if (!is.data.frame(yobs)) yobs <- data.frame(yobs)
      
      # Log likelihood
      ll <- if (is.data.frame(x)) {
        sum(vapply(seq_len(ncol(x)), function(i) {
          do.call(private$.likelihood$log.likelihood, c(list(yobs[[i]], x[[i]]), pars_l))
        }, numeric(1)))
      } else {
        # x is a vector (single series); use first observation column
        yy <- yobs[[1]]
        do.call(private$.likelihood$log.likelihood, c(list(yy, x), pars_l))
      }
      
      ll + lp
    },
    
    #' Run the sampler.
    .calibrate = function(guess, formula, fixed, priors, ...) {
      if (!requireNamespace("mcmc", quietly = TRUE))
        stop("package 'mcmc' is required for Metrop")
      if (!requireNamespace("coda", quietly = TRUE))
        stop("package 'coda' is required for Metrop")
      
      res <- mcmc::metrop(
        obj = private$objective,
        initial = guess,
        ...,
        formula = formula,
        fixed = fixed,
        priors = priors
      )
      colnames(res$batch) <- names(guess)
      res
    },
    
    #' Summarize posterior samples (mean and 95% interval).
    interpret = function(result) {
      s <- summary(coda::as.mcmc(result$batch), quantiles = c(0.025, 0.975))
      as.data.frame(cbind(mean = s$statistics[, "Mean"], s$quantiles))
    }
  ),
  
  public = list(
    #' @description Initialize a Metrop calibrator.
    #'
    #' @param model Model to calibrate.
    #' @param time Time specification (numeric vector or column name in `data`).
    #' @param data Observed data (data.frame).
    #' @param likelihood A `Distribution` object specifying the likelihood.
    #' @param ... Mapping formulas (e.g., `cases ~ I`).
    #' @param cumulative Whether `data` represents cumulative curves.
    #' @param mapping Optional named character vector mapping data columns to model variables.
    initialize = function(model, time, data, likelihood, ..., cumulative = FALSE, mapping = character()) {
      if (!inherits(likelihood, "Distribution"))
        stop("likelihood must be a Distribution object")
      private$.likelihood <- likelihood
      super$initialize(model, time, data, ..., cumulative = cumulative, mapping = mapping)
    }
  ),
  
  active = list(
    #' @field parameters Names of model parameters and likelihood parameters
    parameters = function() c(private$.model$parameters, private$.likelihood$par),
    
    #' @field samples Posterior samples as a `coda::mcmc` object
    samples = function() coda::as.mcmc(private$.details$batch)
  )
)
