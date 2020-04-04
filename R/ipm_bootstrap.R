# Contents of file `ipm_bootstrap.R`

#' Using parametric bootstrap to compute credible interval.
#'
#' @param result an mcmcrs object, the result from ipm_analyse or ipm_analyse_l.
#' @param real an list object, specifying the real value using for simulation.
#' @param priors a string of code to set the prior.
#' @param times a scalar, specifying the number of samples.
#' @param maxtime a scalar, specifying the maximum time to spend on analysis.
#' @param unit a character string specifying the units of time for \code{max.time}. See \code{difftime}.
#' @param save a scalar, the number of (potentially thinned) samples to save.
#' @param chain a scalar, the number of MCMC chains.
#' @return a mcmcrs object, the MCMC outputs.
#' @export


ipm_bootstrap <- function(result,
                          real=NULL,
                          priors=NULL,
                          times=1000,
                          maxtime=10,
                          unit="mins",
                          save=20000L,
                          chain=3){

  ## Simulation
  summ <- sma_summarise(result, measures = "mean")
  p_real <- summ$mcmcr1$mean$p
  f_real <- summ$mcmcr1$mean$f
  N.1_real <- round(summ$mcmcr1$mean$N1)
  sigma_real <- summ$mcmcr1$mean$sigma
  phi_real <- summ$mcmcr1$mean$phi
  xi_real <- summ$mcmcr1$mean$xi
  K <- length(f_real)

  samp <- rep(NA, times)
  j <- 0
  while(j<times){
    flag <- FALSE
    samp[(j+1)] <- tryCatch(ipm_sim_l(K=K, N.1=N.1_real, sigma=sigma_real, p=p_real, xi=xi_real, f=f_real, phi=phi_real),
                            error = function(e){flag <<- TRUE})
    if(flag){
      next
    }else{
      j <- j + 1
    }
  }

  ## Analyse
  res <- rep(NA, times)
  for(i in 1:times){
    res[i] <- ipm_analyse_l(data=nlists(samp[i][[1]]),priors=priors)
  }

  ## Compute the Credible Interval
  real_mat <- c(f_real,N.1_real,p_real,phi_real,sigma_real,xi_real)
  n <- 2*(K-1)+2*K+2
  postmean <- re <- matrix(rep(NA, times*n),ncol=n,nrow=times)
  colnames(postmean) <- colnames(re) <- c(paste("f[",1:K,"]",sep=""),"N1",paste("p[",1:(K-1),"]",sep=""),
                                          paste("phi[",1:(K-1),"]",sep=""),"sigma",paste("xi[",1:K,"]",sep=""))

  for(i in 1:times){
    postmean[i,] <- coef(res[i][[1]])[,2]
    re[i,] <- postmean[i,]-real_mat
  }

  ci.ub <- real_mat - apply(re,2,function(x) quantile(x,0.025))
  ci.lb <- real_mat - apply(re,2,function(x) quantile(x,0.975))

  ## Compute RMSE and Coverage Probability
  rmse <- ipm_rmse(res,times=times)
  coverage <- ipm_coverage(res, times=times)

}


