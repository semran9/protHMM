test_that("the LAG, Bigrams and fusion vectors are all of expected length", {
  expect_equal(length(chmm(system.file("extdata", "1DLHA2-7", package="protHMM"))[[1]]), 800)
  expect_equal(length(chmm(system.file("extdata", "1DLHA2-7", package="protHMM"))[[2]]), 400)
  expect_equal(length(chmm(system.file("extdata", "1DLHA2-7", package="protHMM"))[[3]]), 400)
})
