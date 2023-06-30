test_that("the trigrams function returns the expected vector of length 8000", {
  expect_equal(length(hmm_trigrams(system.file("extdata", "1DLHA2-7", package="protHMM"))), 8000)
})
