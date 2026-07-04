test_that("Model keeps indexed equations and exposes flat equations", {
  groups <- c("A", "K")

  model <- Model$new(
    S[i] ~ -Sum(beta[i, j] * I[j], j = groups) * S[i],
    I[i] ~ Sum(beta[i, j] * I[j], j = groups) * S[i] - gamma[i] * I[i],
    .index = list(i = groups, j = groups)
  )

  expect_true(model$is_stratified())
  expect_equal(model$compartments, c("S_A", "S_K", "I_A", "I_K"))
  expect_equal(model$parameters, c("beta", "gamma"))
  expect_equal(deparse1(model$equations$equations[[1]][[2]]), "S[i]")

  flat <- model$flat_equations()$equations
  expect_named(flat, c("S_A", "S_K", "I_A", "I_K"))

  beta <- matrix(
    c(0.4, 0.2,
      0.3, 0.6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(groups, groups)
  )
  gamma <- c(A = 0.2, K = 0.25)
  env <- list2env(list(beta = beta, gamma = gamma, S_A = 10, I_A = 1, I_K = 2))

  expect_equal(eval(flat$S_A[[3]], envir = env), -8)
  expect_equal(unname(eval(flat$I_A[[3]], envir = env)), 7.8)
})

test_that("symbolic sum and prod expand over model indices", {
  groups <- c("A", "K")

  model <- Model$new(
    X[i] ~ sum(a[j], j = groups) * Prod(p[j], j = groups) * X[i],
    .index = list(i = groups, j = groups)
  )

  flat <- model$flat_equations()$equations
  env <- list2env(list(a = c(A = 2, K = 3), p = c(A = 5, K = 7), X_A = 11))

  expect_equal(unname(eval(flat$X_A[[3]], envir = env)), (2 + 3) * (5 * 7) * 11)
})

test_that("ODE R backend simulates stratified models with structured parameters", {
  skip_if_not_installed("deSolve")

  groups <- c("A", "K")
  model <- Model$new(
    S[i] ~ -Sum(beta[i, j] * I[j], j = groups) * S[i],
    I[i] ~ Sum(beta[i, j] * I[j], j = groups) * S[i] - gamma[i] * I[i],
    .index = list(i = groups, j = groups)
  )

  beta <- matrix(
    c(0.4, 0.2,
      0.3, 0.6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(groups, groups)
  )
  gamma <- c(A = 0.2, K = 0.25)

  out <- ODE$new(model)$simulate(
    0:1,
    y0 = c(S_A = 500, S_K = 500, I_A = 1, I_K = 1),
    parms = list(beta = beta, gamma = gamma)
  )

  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "S_A", "S_K", "I_A", "I_K"))
})

test_that("compiled ODE reports stratification as unsupported for now", {
  skip_if_not_installed("deSolve")

  groups <- c("A", "K")
  model <- Model$new(
    S[i] ~ -Sum(beta[i, j] * I[j], j = groups) * S[i],
    .index = list(i = groups, j = groups)
  )

  expect_error(
    ODE$new(model, compile = TRUE),
    "compiled ODE does not yet support stratified models"
  )
})
