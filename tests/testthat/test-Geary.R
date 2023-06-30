test_that("multiplication works", {
  expect_equal(length(hmm_GA(system.file("extdata", "1DLHA2-7", package="protHMM"))), 180)
})
