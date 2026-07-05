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
  positional <- model$flat_equations(index_mode = "position")$equations
  expect_match(deparse1(positional$S_A[[3]]), "beta_A_A")
  expect_match(deparse1(positional$S_A[[3]]), "beta_A_K")

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

test_that("constructor-level index arguments are consumed as index sets", {
  groups <- c("A", "K")

  model <- Model$new(
    S[i] ~ -Sum(beta[i, j] * I[j], j = groups) * S[i],
    I[i] ~ Sum(beta[i, j] * I[j], j = groups) * S[i] - gamma * I[i],
    i = groups
  )

  expect_true(model$is_stratified())
  expect_equal(model$compartments, c("S_A", "S_K", "I_A", "I_K"))
  expect_equal(model$parameters, c("beta", "gamma"))
  expect_equal(deparse1(model$equations$equations[[2]][[2]]), "I[i]")

  flat <- model$flat_equations()$equations
  expect_named(flat, c("S_A", "S_K", "I_A", "I_K"))
  expect_match(deparse1(flat$I_A[[3]]), "gamma \\* I_A")
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

  ode_body <- deparse1(body(ODE$new(model)$model))
  expect_match(ode_body, "beta_A_A")
  expect_no_match(ode_body, 'beta\\["A", "A"\\]')
})

test_that("positional simulator indexing validates parameter strata order", {
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
  y0 <- c(S_A = 500, S_K = 500, I_A = 1, I_K = 1)

  beta_bad <- beta
  dimnames(beta_bad)[[1]] <- rev(groups)
  expect_error(
    ODE$new(model)$simulate(0:1, y0 = y0, parms = list(beta = beta_bad, gamma = gamma)),
    "parameter beta dimension 1 names must be in model strata order"
  )

  gamma_bad <- gamma[rev(groups)]
  expect_error(
    ODE$new(model)$simulate(0:1, y0 = y0, parms = list(beta = beta, gamma = gamma_bad)),
    "parameter gamma dimension 1 names must be in model strata order"
  )
})

test_that("ODE accepts structured initial values and scalar parameters", {
  skip_if_not_installed("deSolve")

  groups <- c("A", "K")
  model <- Model$new(
    S[i] ~ -Sum(beta[i, j] * I[j], j = groups) * S[i],
    I[i] ~ Sum(beta[i, j] * I[j], j = groups) * S[i] - gamma * I[i],
    i = groups
  )

  sim <- ODE$new(model)
  t <- seq(0, 1, by = 0.1)
  y0 <- list(
    S = c(A = 0.5, K = 0.5),
    I = c(A = 1e-4, K = 1e-4)
  )
  beta <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
  gamma <- 0.2

  out <- sim$simulate(t, y0 = y0, parms = list(beta = beta, gamma = gamma))

  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "S_A", "S_K", "I_A", "I_K"))

  body_text <- deparse1(body(sim$model))
  expect_match(body_text, "gamma \\* I_A")
  expect_match(body_text, "gamma \\* I_K")
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
