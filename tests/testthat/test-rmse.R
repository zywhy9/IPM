context("rmse")

test_that("rmse is correct",{
  set.seed(10)
  data <- ipm_sim()
  result <- ipm_analyse(data,maxtime=1,unit="secs",save=10,chain=1)
  rmse <- ipm_rmse(result,times=1)

  out <- as.matrix(coda::as.mcmc.list(result$mcmcr1))
  real <- c(rep(2,4),500,rep(0.8,3),rep(0.9,3),30,rep(0.8,4))

  for(i in 1:16){
    expect_true(rmse[i]==apply(apply(out,1,function(x){(x-real)^2}),1,function(x){sqrt(mean(x))})[i])
  }
})
