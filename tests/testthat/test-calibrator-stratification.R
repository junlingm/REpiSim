ProbeStratifiedCalibrator <- R6::R6Class(
  "ProbeStratifiedCalibrator",
  inherit = Calibrator,
  private = list(
    simulator = function(model) {
      ODE$new(model)
    },
    .calibrate = function(guess, formula, fixed, ...) {
      list(
        guess = guess,
        formula = formula,
        fixed = fixed,
        simulated = private$simulate(guess, formula, fixed, ...)
      )
    },
    interpret = function(results) {
      results
    }
  )
)

test_that("Calibrator accepts structured stratified initial values, parameters, and guess", {
  skip_if_not_installed("deSolve")

  groups <- c("A", "K")
  model <- Model$new(
    S[i] ~ -Sum(beta[i, j] * I[j], j = groups) * S[i],
    I[i] ~ Sum(beta[i, j] * I[j], j = groups) * S[i] - gamma[i] * I[i],
    .index = list(i = groups, j = groups)
  )

  fit <- ProbeStratifiedCalibrator$new(
    model,
    time = "time",
    data = data.frame(time = c(0, 0.1), obs = c(1e-4, 1.1e-4)),
    obs ~ I_A
  )

  result <- fit$calibrate(
    initial.values = list(
      S = c(A = 0.5, K = 0.5),
      I = c(A = NA, K = 1e-4)
    ),
    parms = list(
      beta = matrix(c(1.5, NA, NA, 2), nrow = 2),
      gamma = c(A = 0.2, K = 0.25)
    ),
    guess = list(
      I = c(A = 1e-4),
      beta = matrix(c(NA, 1.6, 1.7, NA), nrow = 2)
    )
  )

  expect_named(result$guess, c("I_A", "beta_K_A", "beta_A_K"))
  expect_equal(unname(result$guess), c(1e-4, 1.6, 1.7))
  expect_named(
    result$fixed,
    c("S_A", "S_K", "I_K", "beta_A_A", "beta_K_K", "gamma_A", "gamma_K")
  )
  expect_s3_class(result$simulated, "data.frame")
  expect_named(result$simulated, "I_A")
})

test_that("Calibrator accepts flat stratified parameter component names", {
  skip_if_not_installed("deSolve")

  groups <- c("A", "K")
  model <- Model$new(
    S[i] ~ -Sum(beta[i, j] * I[j], j = groups) * S[i],
    I[i] ~ Sum(beta[i, j] * I[j], j = groups) * S[i] - gamma * I[i],
    i = groups
  )

  fit <- ProbeStratifiedCalibrator$new(
    model,
    time = "time",
    data = data.frame(time = c(0, 0.1), obs = c(1e-4, 1.1e-4)),
    obs ~ I_A
  )

  result <- fit$calibrate(
    initial.values = c(S_A = 0.5, S_K = 0.5, I_A = 1e-4, I_K = 1e-4),
    parms = list(beta_A_A = 1.5, beta_K_K = 2, gamma = 0.2),
    guess = c(beta_K_A = 1.6, beta_A_K = 1.7)
  )

  expect_named(result$guess, c("beta_K_A", "beta_A_K"))
  expect_named(result$fixed, c("S_A", "S_K", "I_A", "I_K", "beta_A_A", "beta_K_K", "gamma"))
  expect_s3_class(result$simulated, "data.frame")
})
