context("simulation")

test_that("result is an nlists object",{
  set.seed(10)
  result_cl <- ipm_sim_cl()
  result_l <- ipm_sim_l()
  expect_true(class(result_cl)=="nlists")
  expect_true(class(result_l)=="nlists")
})

test_that("dimensions should be consistent",{
  set.seed(10)
  expect_error(ipm_sim_cl(K=5))
  expect_error(ipm_sim_l(K=3))
})
