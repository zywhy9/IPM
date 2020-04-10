context("summary")

test_that("Summary is correct",{
  set.seed(10)
  data <- ipm_sim_l()
  result <- ipm_analyse(data,maxtime=1,unit="secs",save=10,chain=1)
  sum1 <- ipm_summary(result)

  out <- as.matrix(coda::as.mcmc.list(result$mcmcr1))
  out <- out[,colnames(sum1)]

  for(i in 1:10){
    expect_true(sum1[1,i]==apply(out,2,mean)[i])
    expect_true(sum1[2,i]==apply(out,2,function(x)quantile(x,0.025))[i])
    expect_true(sum1[3,i]==apply(out,2,function(x)quantile(x,0.975))[i])
    expect_true(sum1[4,i]==apply(out,2,var)[i])
  }
})
