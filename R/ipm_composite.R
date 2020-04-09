# Contents of file `ipm_composite.R`

#' Summary for composite likelihood model
#'
#' Estimate the coverage probability and RMSE for composite likelihood model.
#'
#' @param result an mcmcrs object, the result from ipm_analyse.
#' @param data an nlists object, simulation dataset returned by ipm_sim_l.
#' @param real an list object, specifying the real value using for simulation.
#' @param priors a string of code to set the prior.
#' @param times a scalar, specifying the number of samples.
#' @param maxtime a scalar, specifying the maximum time to spend on analysis.
#' @param unit a character string specifying the units of time for \code{max.time}. See \code{difftime}.
#' @param save a scalar, the number of (potentially thinned) samples to save.
#' @param chain a scalar, the number of MCMC chains.
#' @return a list object, containing credible interval, RMSE and coverage probability.
#' @export


ipm_composite <- function(result,
                          data,
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
  K <- length(f_real)
  R <- rep(NA,K-1)
  for(i in 1:(K-1)){
    R[i] <- data[[1]]$Sr[i,i]
  }
  samp <- rep(NA, times)
  j <- 0
  while(j<times){
    flag <- FALSE
    samp[(j+1)] <- tryCatch(ipm_sim_cl(K=K, N.1=N.1_real, sigma=sigma_real, p=p_real, f=f_real, phi=phi_real, R=R),
                            error = function(e){flag <<- TRUE})
    if(flag){
      next
    }else{
      j <- j + 1
    }
  }
  ## Analyse
  res <- rep(NA, times)
  cl <- parallel::makeCluster(spec = parallel::detectCores(logical = FALSE)) # create the parallel cluster
  doParallel::registerDoParallel(cl) # inform doParallel

  res[1:(times/2)] <- foreach(i=1:(times/2),.combine = "c") %dopar% {
    IPM::ipm_analyse_cl(data=nlist::nlists(samp[i][[1]]),
                        priors=priors,maxtime=maxtime,unit=unit,save=save,chain=chain)$result
  }
  res[(times/2+1):times] <- foreach(i=(times/2+1):times,.combine = "c") %dopar% {
    IPM::ipm_analyse_cl(data=nlist::nlists(samp[i][[1]]),
                        priors=priors,maxtime=maxtime,unit=unit,save=save,chain=chain)$result
  }
  parallel::stopCluster(cl)


  ## Compute the Credible Interval
  real_mat <- c(f_real,N.1_real,p_real,phi_real,sigma_real)
  n <- 2*(K-1)+K+2
  postmean <- re <- matrix(NA, ncol=n, nrow=times)
  colnames(postmean) <- colnames(re) <- c(paste0("f[",1:K,"]"),"N1",paste0("p[",1:(K-1),"]"),
                                          paste0("phi[",1:(K-1),"]"),"sigma")

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
                   paste0("phi[",1:(K-1),"]"),"sigma")

  ## Compute RMSE and Coverage Probability
  if(is.null(real)){
    real <- c(rep(2,4),500,rep(0.8,3),rep(0.9,3),30)
  }
  rmse <- ipm_rmse(res,times=times, real=real)

  coverage <- ipm_coverage(res, times=times, real=real, type="cl")
  coverage1 <- ipm_coverage(res, real=real_mat, times=times, type="cl")

  out <- list(var=evar, lower=ci.lb, upper=ci.ub, rmse=rmse, coverage=coverage, coverage1=coverage1)
  return(out)
}


