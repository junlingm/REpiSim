test_that("Compartmental expands stratified transitions into flat events", {
  groups <- c("A", "K")
  model <- Compartmental$new(.index = list(i = groups, j = groups))

  model$transition(
    I[i] <- S[i] ~ Sum(beta[i, j] * I[j], j = groups),
    percapita = TRUE,
    name = "infection"
  )
  model$transition(
    R[i] <- I[i] ~ gamma[i],
    percapita = TRUE,
    name = "recovery"
  )

  expect_equal(model$compartments, c("S_A", "I_A", "S_K", "I_K", "R_A", "R_K"))
  expect_named(model$transitions, c("infection_A", "infection_K", "recovery_A", "recovery_K"))
  expect_equal(names(model$parameter_dimensions()), c("beta", "gamma"))

  infection_a <- model$transitions$infection_A
  expect_equal(infection_a$from, "S_A")
  expect_equal(infection_a$to, "I_A")
  expect_match(deparse1(infection_a$rate), "beta_A_A")
  expect_match(deparse1(infection_a$rate), "S_A")
})

test_that("Compartmental supports transition-local endpoint indices", {
  groups <- c("A", "K")
  model <- Compartmental$new()

  model$transition(
    I[i] <- S[i] ~ Sum(beta[i, j] * I[j], j = groups) * S[i],
    i = groups,
    name = "infection"
  )

  expect_equal(model$compartments, c("S_A", "I_A", "S_K", "I_K"))
  expect_named(model$transitions, c("infection_A", "infection_K"))
  expect_null(model$substitutions$i)
  expect_equal(model$compartment_dimensions()$S, list(groups))
  expect_equal(model$compartment_dimensions()$I, list(groups))
  expect_equal(model$parameter_dimensions()$beta, list(groups, groups))

  infection_a <- model$transitions$infection_A
  expect_equal(infection_a$from, "S_A")
  expect_equal(infection_a$to, "I_A")
  expect_match(deparse1(infection_a$rate), "beta_A_A")
  expect_match(deparse1(infection_a$rate), "S_A")
})

test_that("bare parameters in stratified transitions remain scalar", {
  groups <- c("A", "K")
  model <- Compartmental$new()

  model$transition(R[i] <- I[i] ~ gamma * I[i], i = groups)

  expect_equal(model$parameter_dimensions(), list())
  expect_equal(model$parameters, "gamma")
  expect_equal(model$flat_parameters(), "gamma")
  expect_match(deparse1(model$transitions[["I_A->R_A"]]$rate), "gamma \\* I_A")
})

test_that("Compartmental requires endpoint indices to be specified", {
  groups <- c("A", "K")
  model <- Compartmental$new()

  expect_error(
    model$transition(
      I[i] <- S[i] ~ Sum(beta[i, j] * I[j], j = groups) * S[i],
      name = "infection"
    ),
    "index i must be specified for this transition",
    fixed = TRUE
  )
})

test_that("ODE simulates stratified Compartmental models with structured inputs", {
  skip_if_not_installed("deSolve")

  groups <- c("A", "K")
  model <- Compartmental$new(.index = list(i = groups, j = groups))
  model$transition(I[i] <- S[i] ~ Sum(beta[i, j] * I[j], j = groups), percapita = TRUE)
  model$transition(R[i] <- I[i] ~ gamma[i], percapita = TRUE)

  out <- ODE$new(model)$simulate(
    0:1,
    y0 = list(
      S = c(A = 500, K = 500),
      I = c(A = 1, K = 1),
      R = c(A = 0, K = 0)
    ),
    parms = list(
      beta = matrix(
        c(0.4, 0.2,
          0.3, 0.6),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(groups, groups)
      ),
      gamma = c(A = 0.2, K = 0.25)
    )
  )

  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "S_A", "S_K", "I_A", "I_K", "R_A", "R_K"))
})

test_that("Gillespie sees stratified Compartmental transitions as separate events", {
  groups <- c("A", "K")
  model <- Compartmental$new(.index = list(i = groups))
  model$transition(NULL <- I[i] ~ gamma[i], percapita = TRUE, name = "death")

  expect_equal(model$compartments, c("I_A", "I_K"))
  expect_named(model$transitions, c("death_A", "death_K"))

  sim <- Gillespie$new(model, compile = FALSE)
  out <- sim$simulate(
    0:3,
    y0 = list(I = c(A = 2, K = 2)),
    parms = list(gamma = c(A = 100, K = 100))
  )

  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "I_A", "I_K"))
})
