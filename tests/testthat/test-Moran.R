test_that("the Moran autocorrelation vector is of correct length 180", {
  expect_equal(length(hmm_MA(system.file("extdata", "1DLHA2-7", package="protHMM"))), 180)
})
