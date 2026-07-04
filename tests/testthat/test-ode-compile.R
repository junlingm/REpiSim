test_that("ODE compile=TRUE uses a compiled model compatible with deSolve", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("cpp11")

  model <- Model$new(S ~ -exp(a) * S)
  y0 <- c(S = 10)
  parms <- c(a = log(0.1))

  r_model <- ODE$new(model, compile = FALSE)
  c_model <- ODE$new(model, compile = TRUE)

  r_out <- r_model$simulate(0:3, y0 = y0, parms = parms)
  c_out <- c_model$simulate(0:3, y0 = y0, parms = parms)

  expect_equal(c_out, r_out, tolerance = 1e-7, ignore_attr = TRUE)
  expect_type(c_model$compiled, "list")
  expect_match(c_model$compiled$rhs, "cpp11::doubles ode_rhs", fixed = TRUE)

  value <- c_model$model(0, y0, parms)
  expect_type(value, "list")
  expect_length(value, 2)
  expect_length(value[[1]], 1)
})

test_that("compiled ODE models return substitutions like the R backend", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("cpp11")

  model <- Compartmental$new(S, I, R)
  model$transition(S -> I ~ beta * S * I / N, N = S + I + R, name = "infection")
  model$transition(I -> R ~ gamma * I, name = "recovery")

  y0 <- c(S = 1000, I = 10, R = 0)
  parms <- c(beta = 0.4, gamma = 0.2)
  vars <- c("time", "S", "I", "R", "N")

  r_out <- ODE$new(model)$simulate(0:5, y0 = y0, parms = parms, vars = vars)
  c_out <- ODE$new(model, compile = TRUE)$simulate(0:5, y0 = y0, parms = parms, vars = vars)

  expect_equal(c_out, r_out, tolerance = 1e-7, ignore_attr = TRUE)
})

test_that("compiled ODE models handle common scalar math helpers", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("cpp11")

  model <- Model$new(S ~ -ifelse(S > 0, max(0, exp(a)) * sign(S), 0))
  y0 <- c(S = 10)
  parms <- c(a = log(0.1))

  r_out <- ODE$new(model)$simulate(0:3, y0 = y0, parms = parms)
  c_out <- ODE$new(model, compile = TRUE)$simulate(0:3, y0 = y0, parms = parms)

  expect_equal(c_out, r_out, tolerance = 1e-7, ignore_attr = TRUE)
})
