test_that("hdsbm_fit converges and produces correct dimensions", {
  data(toy_hdsbm)
  fit <- hdsbm_fit(toy_hdsbm$Y, k = 3, nrep = 1, maxit = 10, verbose = FALSE)
  expect_s3_class(fit, "hdsbm")
  expect_equal(dim(fit$best_fit$Psi), c(3, 3))
  expect_equal(dim(fit$best_fit$Pi), c(3, 3))
})