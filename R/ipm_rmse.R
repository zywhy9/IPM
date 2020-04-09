# Contents of file `ipm_rmse.R`

#' Root mean square error
#'
#' Calculate the root mean square error.
#'
#' @param result an mcmcrs object, the result from ipm_analyse or ipm_analyse_l.
#' @param real an list object, specifying the real value using for simulation.
#' @param times a scalar, specifying the number of samples.
#' @return a vector, the RMSE of each parameter.
#' @export


ipm_rmse <- function(result,
                     real=NULL,
                     times=1000){

  ## Compute RMSE
  if(is.null(real)){
    real <- c(rep(2,4),500,rep(0.8,3),rep(0.9,3),30,rep(0.8,4))
  }

  post <- c()
  for(i in 1:times){
    post <- rbind(as.matrix(coda::as.mcmc.list(result[i][[1]])))
  }
  se <- apply(post,1,function(x){(x-real)^2})
  rmse <- apply(se,1,function(x){sqrt(mean(x))})

  return(rmse)
}


