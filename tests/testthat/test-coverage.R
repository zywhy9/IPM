context("coverage")

test_that("coverage is correct",{
  set.seed(10)
  data <- ipm_sim_l()
  result <- ipm_analyse(data,maxtime=1,unit="secs",save=10, chain=1)
  cp <- ipm_coverage(result,times=1)

  out <- as.matrix(coda::as.mcmc.list(result$mcmcr1))
  real <- c(rep(2,4),500,rep(0.8,3),rep(0.9,3),30,rep(0.8,4))
  lb <- apply(out,2,function(x)quantile(x,0.025))
  ub <- apply(out,2,function(x)quantile(x,0.975))

  for(i in 1:16){
    expect_true(cp[i]==as.numeric((real>=lb)&(real<=ub))[i])
  }
})
