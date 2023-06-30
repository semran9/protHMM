test_that("the GSD vector is of correct length", {
  expect_equal(length(hmm_GSD(system.file("extdata", "1DLHA2-7", package="protHMM"))), 300)
})
