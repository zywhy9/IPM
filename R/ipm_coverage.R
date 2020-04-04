# Contents of file `ipm_coverage.R`

#' Coverage probability
#'
#' Calculate the coverage probability.
#'
#' @param result an mcmcrs object, the result from ipm_analyse or ipm_analyse_l.
#' @param real an list object, specifying the real value using for simulation.
#' @param times a scalar, specifying the number of samples.
#' @return a mcmcrs object, the MCMC outputs.
#' @export


ipm_coverage <- function(result,
                     real=NULL,
                     times=1000){

  ## Compute RMSE
  if(is.null(real)){
    real <- c(rep(2,4),500,rep(0.8,3),rep(0.9,3),30,rep(0.8,4))
  }

  K <- dim(result[1][[1]]$f)[3]
  n <- 2*(K-1)+2*K+2
  count <- rep(0,n)
  for(i in 1:times){
    count <- count + as.numeric((real>=coef(result[i][[1]])[,"lower"])&(real<=coef(result[i][[1]])[,"upper"]))
  }

  out <- count/times
  names(out) <- c(paste("f[",1:K,"]",sep=""),"N1",paste("p[",1:(K-1),"]",sep=""),
                  paste("phi[",1:(K-1),"]",sep=""),"sigma",paste("xi[",1:K,"]",sep=""))
  return(out)
}


