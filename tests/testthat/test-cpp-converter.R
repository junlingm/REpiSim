test_that("CppConverter translates arithmetic and indexing expressions", {
  cpp <- REpiSim:::CppConverter$new(
    compartments = c("S", "I"),
    parameters = c("beta"),
    alias = list()
  )

  expect_equal(cpp$expr(quote(beta * S^2 / (S + I))), "beta * pow(S, 2) / (S + I)")
  expect_equal(cpp$expr(quote(exp(beta) * sqrt(S))), "exp(beta) * sqrt(S)")
  expect_equal(cpp$expr(quote(abs(-S) + round(beta))), "fabs(-S) + round(beta)")
  expect_equal(cpp$expr(quote(S > 0 & I <= 1)), "S > 0 && I <= 1")
  expect_equal(cpp$expr(quote(x[[i]])), "x[i]")
  expect_equal(cpp$expr(quote(x[[i]][[j]])), "x[i][j]")
})

test_that("CppConverter translates common scalar math helpers", {
  cpp <- REpiSim:::CppConverter$new(
    compartments = c("S", "I"),
    parameters = c("beta"),
    alias = list()
  )

  expect_equal(cpp$expr(quote(sign(S))), "((S > 0) - (S < 0))")
  expect_equal(cpp$expr(quote(min(S, I, beta))), "fmin(fmin(S, I), beta)")
  expect_equal(cpp$expr(quote(max(S, I, beta))), "fmax(fmax(S, I), beta)")
  expect_equal(cpp$expr(quote(sum(S, I, beta))), "(S + I + beta)")
  expect_equal(cpp$expr(quote(prod(S, I, beta))), "(S * I * beta)")
  expect_equal(cpp$expr(quote(all(is.finite(S), I > 0))), "(std::isfinite(S) && I > 0)")
  expect_equal(cpp$expr(quote(any(is.nan(S), is.infinite(I)))), "(std::isnan(S) || std::isinf(I))")
  expect_equal(cpp$expr(quote(ifelse(S > 0, S, 0))), "((S > 0) ? (S) : (0))")
})

test_that("CppConverter errors clearly for non-scalar whitelisted helpers", {
  cpp <- REpiSim:::CppConverter$new(
    compartments = c("S", "I"),
    parameters = c("beta"),
    alias = list()
  )

  expect_error(
    cpp$expr(quote(mean(S))),
    "allowed in R models but is not supported by the compiled C\\+\\+ converter"
  )
})

test_that("CppConverter emits model context bindings in order", {
  cpp <- REpiSim:::CppConverter$new(
    compartments = c("S", "I"),
    parameters = c("beta"),
    alias = list(N = quote(N <- S + I), incidence = quote(incidence <- beta * S * I / N))
  )

  expect_equal(
    cpp$bind_vars(c("S", "I"), "y"),
    list("double S = __y[0];", "double I = __y[1];")
  )
  expect_equal(
    cpp$substitutions(),
    list(N = "double N = S + I;", incidence = "double incidence = beta * S * I / N;")
  )
})
