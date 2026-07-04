test_that("CppConverter translates arithmetic and indexing expressions", {
  cpp <- REpiSim:::CppConverter$new(
    compartments = c("S", "I"),
    parameters = c("beta"),
    alias = list()
  )

  expect_equal(cpp$expr(quote(beta * S^2 / (S + I))), "beta * pow(S, 2) / (S + I)")
  expect_equal(cpp$expr(quote(x[[i]])), "x[i]")
  expect_equal(cpp$expr(quote(x[[i]][[j]])), "x[i][j]")
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
