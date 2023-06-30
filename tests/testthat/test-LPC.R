test_that("the LPC vector is of expected length 280", {
  expect_equal(length(hmm_LPC(system.file("extdata", "1DLHA2-7", package="protHMM"))), 280)
})
