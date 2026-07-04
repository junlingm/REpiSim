ProbeCalibrator <- R6::R6Class(
  "ProbeCalibrator",
  inherit = Calibrator,

  private = list(
    simulator = function(model) {
      ODE$new(model, compile = private$.compile)
    }
  ),

  public = list(
    simulator_is_compiled = function() {
      !is.null(private$.simulator$compiled)
    }
  )
)

test_that("Calibrator forwards compile to default simulators", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("cpp11")

  model <- Model$new(S ~ -a * S)
  data <- data.frame(time = 0:2, obs = c(10, 9, 8))

  r_fit <- ProbeCalibrator$new(model, time = "time", data = data, obs ~ S)
  c_fit <- ProbeCalibrator$new(model, time = "time", data = data, obs ~ S, compile = TRUE)

  expect_false(r_fit$simulator_is_compiled())
  expect_true(c_fit$simulator_is_compiled())
})

test_that("LeastSquare accepts compile=TRUE", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("cpp11")

  model <- Model$new(S ~ -a * S)
  data <- data.frame(time = 0:2, obs = c(10, 9, 8))

  fit <- LeastSquare$new(model, time = "time", data = data, obs ~ S, compile = TRUE)

  expect_equal(fit$parameters, "a")
})
