test_that("icl function works on hdsbm objects and raw outputs", {
  data(toy_hdsbm)
  fit <- hdsbm_fit(toy_hdsbm$Y, k = 3, nrep = 1, maxit = 10, verbose = FALSE)
  
  # Test single numeric return
  val <- icl(fit)
  expect_type(val, "double")
  expect_equal(val, fit$icl$ICL)
  
  # Test detailed return
  det <- icl(fit, detailed = TRUE)
  expect_type(det, "list")
  expect_named(det, c("ICL", "lk", "pen_total", "pen_components"))
  
  # Test backward-compatible icl_dyn_hyb
  compat <- icl_dyn_hyb(fit$best_fit)
  expect_equal(compat$ICL, val)
})