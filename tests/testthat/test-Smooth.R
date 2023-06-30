test_that("The matrix returned by the Smoothing function is of proper dimensions", {
  expect_equal(nrow(hmm_smooth(system.file("extdata", "1DLHA2-7", package="protHMM"))),
               nrow(hmm_read(system.file("extdata", "1DLHA2-7", package="protHMM"))))
  expect_equal(ncol(hmm_smooth(system.file("extdata", "1DLHA2-7", package="protHMM"))),
               ncol(hmm_read(system.file("extdata", "1DLHA2-7", package="protHMM"))))
})
