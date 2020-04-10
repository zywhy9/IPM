# Contents of file `ipm_bootstrap.R`

#' Parametric Bootstrap Estimation
#'
#' Using parametric bootstrap to compute credible interval.
#'
#' @param result an mcmcrs object, the result from ipm_analyse or ipm_analyse_l.
#' @param type a string, specifying the model used to analyse.
#' @param real an list object, specifying the real value using for simulation.
#' @param priors a string of code to set the prior.
#' @param times a scalar, specifying the number of samples.
#' @param maxtime a scalar, specifying the maximum time to spend on analysis.
#' @param unit a character string specifying the units of time for \code{max.time}. See \code{difftime}.
#' @param save a scalar, the number of (potentially thinned) samples to save.
#' @param chain a scalar, the number of MCMC chains.
#' @return a list object, containing credible interval, RMSE and coverage probability.
#' @export


ipm_bootstrap <- function(result,
                          type="l",
                          real=NULL,
                          priors=NULL,
                          times=1000,
                          maxtime=5,
                          unit="mins",
                          save=20000L,
                          chain=3){

  ## Simulation
  summ <- simanalyse::sma_summarise(result, measures = "mean")
  p_real <- summ$mcmcr1$mean$p
  f_real <- summ$mcmcr1$mean$f
  N.1_real <- round(summ$mcmcr1$mean$N1)
  sigma_real <- summ$mcmcr1$mean$sigma
  phi_real <- summ$mcmcr1$mean$phi
  xi_real <- summ$mcmcr1$mean$xi
  K <- length(f_real)

  samp <- rep(NA, times)
  j <- 0
  err <- 0
  while(j<times){
    if(err > 10000){
      stop("use simanalyse::sma_summarise(result, measures = \"mean\") to check if any probability is larger than 1.
          If so, try to run ipm_analyse or ipm_analyse_l again.devt")
    }
    flag <- FALSE
    samp[(j+1)] <- tryCatch(ipm_sim_l(K=K, N.1=N.1_real, sigma=sigma_real, p=p_real, xi=xi_real, f=f_real, phi=phi_real),
                            error = function(e){flag <<- TRUE})
    if(flag){
      err <- err + 1
      next
    }else{
      j <- j + 1
    }
  }

  ## Analyse
  res <- rep(NA, times)
  cl <- parallel::makeCluster(spec = parallel::detectCores(logical = FALSE)) # create the parallel cluster
  doParallel::registerDoParallel(cl) # inform doParallel

  if(type=="l"){
    res[1:(times/2)] <- foreach(i=1:(times/2),.combine = "c") %dopar% {
      IPM::ipm_analyse_l(data=nlist::nlists(samp[i][[1]]),
                         priors=priors,maxtime=maxtime,unit=unit,save=save,chain=chain)
    }
    res[(times/2+1):times] <- foreach(i=(times/2+1):times,.combine = "c") %dopar% {
      IPM::ipm_analyse_l(data=nlist::nlists(samp[i][[1]]),
                         priors=priors,maxtime=maxtime,unit=unit,save=save,chain=chain)
    }
  }else{
    res[1:(times/2)] <- foreach(i=1:(times/2),.combine = "c") %dopar% {
      IPM::ipm_analyse(data=nlist::nlists(samp[i][[1]]),
                         priors=priors,maxtime=maxtime,unit=unit,save=save,chain=chain)
    }
    res[(times/2+1):times] <- foreach(i=(times/2+1):times,.combine = "c") %dopar% {
      IPM::ipm_analyse(data=nlist::nlists(samp[i][[1]]),
                         priors=priors,maxtime=maxtime,unit=unit,save=save,chain=chain)
    }
  }



  parallel::stopCluster(cl)


  ## Compute the Credible Interval
  real_mat <- c(f_real,N.1_real,p_real,phi_real,sigma_real,xi_real)
  n <- 2*(K-1)+2*K+2
  postmean <- re <- matrix(NA, ncol=n, nrow=times)
  colnames(postmean) <- colnames(re) <- c(paste0("f[",1:K,"]"),"N1",paste0("p[",1:(K-1),"]"),
                                          paste0("phi[",1:(K-1),"]"),"sigma",paste0("xi[",1:K,"]"))

  for(i in 1:times){
    postmean[i,] <- coef(res[i][[1]])[,2]
    re[i,] <- postmean[i,]-real_mat
  }

  ci.ub <- real_mat - apply(re,2,function(x) quantile(x,0.025))
  ci.lb <- real_mat - apply(re,2,function(x) quantile(x,0.975))

  ## Compute Variance
  postvar <- 0
  for(i in 1:times){
    postvar <- postvar + (coef(res[i][[1]])[,"sd"])^2
  }
  evar <- postvar/times
  names(evar) <- c(paste0("f[",1:K,"]"),"N1",paste0("p[",1:(K-1),"]"),
                   paste0("phi[",1:(K-1),"]"),"sigma",paste0("xi[",1:K,"]"))

  ## Compute RMSE and Coverage Probability
  rmse <- ipm_rmse(res,real=real, times=times)
  coverage <- ipm_coverage(res, real=real, times=times)
  coverage1 <- ipm_coverage(res, real=real_mat, times=times)

  out <- list(var=evar, lower=ci.lb, upper=ci.ub, rmse=rmse, coverage=coverage, coverage1=coverage1)
  return(out)
}


