test_that("Compartmental transitions declare compartments without constructor args", {
  model <- Compartmental$new()

  model$transition(I <- S ~ beta * S * I / N, N = S + I + R, name = "infection")
  model$transition(R <- I ~ gamma * I, name = "recovery")

  expect_equal(model$compartments, c("S", "I", "R"))
  expect_named(model$transitions, c("infection", "recovery"))
  expect_named(model$flat_equations()$equations, c("S", "I", "R"))

  out <- ODE$new(model)$simulate(
    0:1,
    y0 = c(S = 1000, I = 10, R = 0),
    parms = c(beta = 0.4, gamma = 0.2),
    vars = c("S", "I", "R")
  )
  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "S", "I", "R"))
})

test_that("explicit Compartmental declarations still control compartment order", {
  model <- Compartmental$new(I, S, R)
  model$transition(I <- S ~ beta * S * I / N, N = S + I + R)
  model$transition(R <- I ~ gamma * I)

  expect_equal(model$compartments, c("I", "S", "R"))
})

test_that("Compartmental NULL endpoints do not create NULL compartments", {
  model <- Compartmental$new()

  model$transition(I <- NULL ~ birth)
  model$transition(NULL <- I ~ death, percapita = TRUE)

  expect_equal(model$compartments, "I")
  expect_equal(length(model$transitions), 2)
})
