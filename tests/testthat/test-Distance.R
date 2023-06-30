test_that("multiplication works", {
  expect_equal(length(hmm_distance(system.file("extdata", "1DLHA2-7", package="protHMM"), system.file("extdata", "1TEN-7", package="protHMM"))), 1)
  expect_equal(is.double(hmm_distance(system.file("extdata", "1DLHA2-7", package="protHMM"), system.file("extdata", "1TEN-7", package="protHMM"))), TRUE)
})
