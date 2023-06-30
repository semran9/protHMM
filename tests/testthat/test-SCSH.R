test_that("the SCSH2 vectors are both of the expected lengths", {
  expect_equal(length(hmm_SCSH(system.file("extdata", "1DLHA2-7", package="protHMM"))[[1]]), 400)
  expect_equal(length(hmm_SCSH(system.file("extdata", "1DLHA2-7", package="protHMM"))[[2]]), 8000)
})
