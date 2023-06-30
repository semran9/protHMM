test_that("the Separated Dimers vector is of expected length 400", {
  expect_equal(length(hmm_SepDim(system.file("extdata", "1DLHA2-7", package="protHMM"))), 400)
})
