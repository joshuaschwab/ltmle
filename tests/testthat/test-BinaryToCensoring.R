context("Test BinaryToCensoring") 

test_that("BinaryToCensoring works", {
  expect_equal(BinaryToCensoring(is.censored = c(0, 1, 1, 0, NA)), BinaryToCensoring(is.uncensored=c(1, 0, 0, 1, NA)))
}
)