# ==============================================================================
# CGillespie compatibility wrapper
# ==============================================================================
#
# The Gillespie simulator now supports both R and compiled cpp11 backends through
# `Gillespie$new(model, compile = FALSE/TRUE)`. This class is retained so
# existing code that calls `CGillespie$new(model)` keeps working.
# ==============================================================================

#' R6 compatibility class for the compiled Gillespie backend
#'
#' `CGillespie$new(model)` is equivalent to `Gillespie$new(model, compile = TRUE)`.
#'
#' @docType class
#' @export
CGillespie <- R6Class(
  "CGillespie",
  inherit = Gillespie,
  
  public = list(
    #' @description
    #' Construct a compiled Gillespie simulator.
    #'
    #' @param model A [Compartmental] model.
    initialize = function(model) {
      super$initialize(model, compile = TRUE)
    }
  )
)
