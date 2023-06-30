test_that("the bigrams vector is the expected length 400 output", {
  expect_equal(length(hmm_bigrams(system.file("extdata", "1DLHA2-7", package="protHMM"))), 400)
})
