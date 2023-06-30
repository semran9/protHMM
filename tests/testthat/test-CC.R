test_that("the CC vector is of expected length, with both default and changed lag values", {
  expect_equal(length(hmm_cc(system.file("extdata", "1DLHA2-7", package="protHMM"))), 1520)
  expect_equal(length(hmm_cc(system.file("extdata", "1DLHA2-7", package="protHMM"), 5)), 1900)
})
