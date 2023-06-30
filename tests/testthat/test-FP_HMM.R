test_that("the lengths of the two vectors given by the fp_hmm function are the expected lengths", {
  expect_equal(length(fp_hmm(system.file("extdata", "1DLHA2-7", package="protHMM"))[[1]]), 20)
  expect_equal(length(fp_hmm(system.file("extdata", "1DLHA2-7", package="protHMM"))[[2]]), 400)
})
