context("analysis")

test_that("result is a mcmcrs object",{
  set.seed(10)
  data <- ipm_sim_l()
  result_cl <- ipm_analyse(data,maxtime=1,unit="secs",save=10)
  result_l <- ipm_analyse_l(data,maxtime=1,unit="secs",save=10)
  expect_true(class(result_cl)=="mcmcrs")
  expect_true(class(result_l)=="mcmcrs")
})
