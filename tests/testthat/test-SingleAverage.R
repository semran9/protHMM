test_that("the Single Average vector is returned with the correct length of 400", {
  expect_equal(length(hmm_Single_Average(system.file("extdata", "1DLHA2-7", package="protHMM"))), 400)
})
