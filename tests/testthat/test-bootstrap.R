context("bootstrap")

test_that("result is a list object",{
  set.seed(10)
  data <- ipm_sim_l()
  result <- ipm_analyse(data,maxtime=1,unit="secs",save=10)
  out_l <- ipm_bootstrap(result,maxtime=1,unit="secs",save=10,times=10)
  out_cl <- ipm_composite(result,data,maxtime=1,unit="secs",save=10,times=10)
  expect_true(class(out_cl)=="list")
  expect_true(class(out_l)=="list")
})
