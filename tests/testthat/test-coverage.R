context("coverage")

test_that("result is a numeric object",{
  set.seed(10)
  data <- ipm_sim_l()
  result <- ipm_analyse(data,maxtime=1,unit="secs",save=10)
  out <- ipm_coverage(result,times=1)
  expect_true(class(out)=="numeric")
})
