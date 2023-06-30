test_that("multiplication works", {
  expect_equal(length(hmm_LBP(system.file("extdata", "1DLHA2-7", package="protHMM"))), 256)
})
