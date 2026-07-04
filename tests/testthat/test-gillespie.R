pure_death_model <- function() {
  model <- Compartmental$new(I)
  model$transition(NULL <- I ~ gamma, percapita = TRUE, name = "death")
  model
}

test_that("Gillespie R backend simulates a pure-death process to extinction", {
  model <- pure_death_model()
  sim <- Gillespie$new(model, compile = FALSE)

  set.seed(1)
  out <- sim$simulate(0:3, y0 = c(I = 3), parms = c(gamma = 100))

  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "I"))
  expect_equal(out$time, 0:3)
  expect_equal(out$I, c(3, 0, 0, 0))
})

test_that("RGillespie remains a compile=FALSE compatibility wrapper", {
  model <- pure_death_model()
  sim <- RGillespie$new(model)

  expect_null(sim$compiled)
})

test_that("Gillespie compiled backend matches the pure-death smoke case", {
  skip_if_not_installed("cpp11")

  model <- pure_death_model()
  sim <- Gillespie$new(model, compile = TRUE)

  set.seed(1)
  out <- sim$simulate(0:3, y0 = c(I = 3), parms = c(gamma = 100))

  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "I"))
  expect_equal(out$time, 0:3)
  expect_equal(out$I, c(3, 0, 0, 0))
  expect_match(sim$compiled$program, "inline cpp11::doubles step", fixed = TRUE)
})

test_that("CGillespie remains a compile=TRUE compatibility wrapper", {
  skip_if_not_installed("cpp11")

  model <- pure_death_model()
  sim <- CGillespie$new(model)

  expect_type(sim$compiled, "list")
  expect_true(is.function(sim$compiled$simulate))
})
