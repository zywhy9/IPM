context("analysis/bootstrap")

test_that("analysis is correct",{
  set.seed(10)
  ## Simulate data
  n_obs <- 5
  mu <- rnorm(1)
  X <- rnorm(n_obs,mean=mu)

  # simulate some parameter values
  n_sim <- 15
  Pars <- replicate(n = n_sim, expr = {
    list(mu = rnorm(1))
  }, simplify = FALSE)

  # log-posterior calculation in R
  logpost <- function(mu,X) {
    mu_p <- sum(X)/(1+n_obs)
    sd_p <- sqrt(1/(1+n_obs))
    dnorm(mu,mean=mu_p,sd=sd_p,log = T)
  }

  lp_r <- sapply(1:n_sim, function(ii) {
    do.call(logpost, c(Pars[[ii]], list(X=X)))
  })

  # log-posterior calculation in IPM
  model <- "
  for(i in 1:5){
   x[i] ~ dnorm(mu,1)
  }
  mu ~ dnorm(0,1)
  mu_p <- sum(x)/6
  logp <- log(1/sqrt(2*pi/6)*exp(-(mu-mu_p)^2*3))
  "
  temp <- c()
  for(ii in 1:n_sim){
    data_a <- nlists(nlist(x=X,mu=Pars[[ii]]$mu,pi=pi))
    result <- ipm_analyse(data=data_a, model=model,priors="",moni="logp",inits=list(),chain=1,save=10)
    temp <- as.numeric(c(temp,as.matrix(coda::as.mcmc.list(result$mcmcr1$logp))[1,]))
  }

  expect_equal(lp_r, temp, tolerance = 0.000001)
})
