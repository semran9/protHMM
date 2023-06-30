test_that("mthe length of the improed psuedo hmm is as expected", {
  expect_equal(length(IM_psehmm(system.file("extdata", "1DLHA2-7", package="protHMM"))), 189)
})
