test_that("the pseudo hmm vector is of the expected length", {
  expect_equal(length(pse_hmm(system.file("extdata", "1DLHA2-7", package="protHMM"))), 320)
})
