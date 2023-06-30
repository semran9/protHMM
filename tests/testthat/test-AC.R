test_that("the AC vector is the expected length, with both default and changed lag values", {
  expect_equal(length(hmm_ac(system.file("extdata", "1DLHA2-7", package="protHMM"))), 80)
  expect_equal(length(hmm_ac(system.file("extdata", "1DLHA2-7", package="protHMM"), 5)), 100)
})
