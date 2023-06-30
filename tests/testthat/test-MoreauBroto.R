test_that("the Moreau-Broto autocorrelation vector is of correct length 180", {
  expect_equal(length(hmm_MB(system.file("extdata", "1DLHA2-7", package="protHMM"))), 180)
})
