# Contents of file `ipm_summary.R`

#' Summary of the result from analysis
#'
#' @param result an mcmcrs object, the result from ipm_analyse or ipm_analyse_l.
#' @return a matirx, with summary staitstics.
#' @export


ipm_summary <- function(result){

  summ <- sma_summarise(result)

  p.mean <- summ$mcmcr1$mean$p
  f.mean <- summ$mcmcr1$mean$f
  N1.mean <- summ$mcmcr1$mean$N1
  sig.mean <- summ$mcmcr1$mean$sigma
  phi.mean <- summ$mcmcr1$mean$phi

  p.lb <- summ$mcmcr1$lower.q$p
  f.lb <- summ$mcmcr1$lower.q$f
  N1.lb <- summ$mcmcr1$lower.q$N1
  sig.lb <- summ$mcmcr1$lower.q$sigma
  phi.lb <- summ$mcmcr1$lower.q$phi

  p.ub <- summ$mcmcr1$upper.q$p
  f.ub <- summ$mcmcr1$upper.q$f
  N1.ub <- summ$mcmcr1$upper.q$N1
  sig.ub <- summ$mcmcr1$upper.q$sigma
  phi.ub <- summ$mcmcr1$upper.q$phi

  res <- rbind(c(p.mean,f.mean,N1.mean,sig.mean,phi.mean),
               c(p.lb,f.lb,N1.lb,sig.lb,phi.lb),
               c(p.ub,f.ub,N1.ub,sig.ub,phi.ub))

  K <- length(f.mean)

  colnames(res) <- c(paste("p[",1:(K-1),"]",sep=""), paste("f[",1:K,"]",sep=""), "N1", "sigma",
                     paste("phi[",1:(K-1),"]",sep=""))

  rownames(res) <- c("mean","2.5%","97.5%")

  return(res)
}



