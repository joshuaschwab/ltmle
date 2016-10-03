context("General Utilities")

test_that("safe.solve works", {
  a <- matrix(c(1, 2, 2, 4), 2, 2)
  expect_warning(expect_equal(safe.solve(a), matrix(nrow=2, ncol=2)))
  expect_warning(expect_equal(safe.solve(a, matrix(5:6, ncol=1)), matrix(nrow=2, ncol=1)))
})