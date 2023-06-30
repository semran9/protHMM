test_that("the first twenty columns of thehidden markov model file are read in", {
  expect_equal(ncol(hmm_read(system.file("extdata", "1DLHA2-7", package="protHMM"))), 20)
  expect_equal(is.numeric(hmm_read(system.file("extdata", "1DLHA2-7", package="protHMM"))), TRUE)
})
