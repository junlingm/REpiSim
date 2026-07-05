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

test_that("LeastSquare compile=TRUE supports stratified Compartmental parameters", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("cpp11")

  groups <- c("A", "K")
  model <- Compartmental$new()
  model$transition(I[i] <- S[i] ~ Sum(beta[i, j] * I[j], j = groups) * S[i], i = groups)
  model$transition(R[i] <- I[i] ~ gamma * I[i], i = groups)

  t <- seq(0, 1, by = 0.1)
  y0 <- list(
    S = c(A = 0.5, K = 0.5),
    I = c(A = 1e-4, K = 1e-4),
    R = c(A = 0, K = 0)
  )
  beta <- matrix(c(1.5, 0.7, 0.7, 2), nrow = 2)
  data <- ODE$new(model)$simulate(t, y0 = y0, parms = list(beta = beta, gamma = 0.2))

  fit <- LeastSquare$new(
    model,
    data = data,
    time = "time",
    mapping = c(I_A = "I_A", I_K = "I_K"),
    compile = TRUE
  )

  expect_equal(fit$parameters, c("beta_A_A", "beta_A_K", "gamma", "beta_K_A", "beta_K_K"))
})
