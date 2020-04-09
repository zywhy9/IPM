# Contents of file `ipm_coverage.R`

#' Coverage probability
#'
#' Calculate the coverage probability.
#'
#' @param result an mcmcrs object, the result from ipm_analyse or ipm_analyse_l or ipm_analyse_cl.
#' @param real an list object, specifying the real value using for simulation.
#' @param times a scalar, specifying the number of samples.
#' @param type a string, specifying if it's from ipm_analyse_cl.
#' @return a vector, the coverage probability of each parameter.
#' @export


ipm_coverage <- function(result,
                     real=NULL,
                     times=1000,
                     type="l"){

  ## Compute Coverage Probability
  if(is.null(real)){
    real <- c(rep(2,4),500,rep(0.8,3),rep(0.9,3),30,rep(0.8,4))
  }
  if(type=="l"){
    K <- dim(result[1][[1]]$f)[3]
    n <- 2*(K-1)+2*K+2
    count <- rep(0,n)
    for(i in 1:times){
      count <- count + as.numeric((real>=coef(result[i][[1]])[,"lower"])&(real<=coef(result[i][[1]])[,"upper"]))
    }

    out <- count/times
    names(out) <- c(paste0("f[",1:K,"]"),"N1",paste0("p[",1:(K-1),"]"),
                    paste0("phi[",1:(K-1),"]"),"sigma",paste0("xi[",1:K,"]"))
  }else{
    K <- dim(result[1][[1]]$f)[3]
    n <- 2*(K-1)+K+2
    count <- rep(0,n)
    for(i in 1:times){
      count <- count + as.numeric((real>=coef(result[i][[1]])[,"lower"])&(real<=coef(result[i][[1]])[,"upper"]))
    }

    out <- count/times
    names(out) <- c(paste0("f[",1:K,"]"),"N1",paste0("p[",1:(K-1),"]"),
                    paste0("phi[",1:(K-1),"]"),"sigma")
  }

  return(out)
}


