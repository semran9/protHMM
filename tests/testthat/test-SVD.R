test_that("the SVD vector is of expected length 20", {
  expect_equal(length(hmm_svd(system.file("extdata", "1DLHA2-7", package="protHMM"))), 20)
})
