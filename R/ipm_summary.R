# Contents of file `ipm_summary.R`

#' Summary of the result from analysis
#'
#' @param result an mcmcrs object, the result from ipm_analyse or ipm_analyse_l.
#' @return a matirx, with summary staitstics.
#' @export


ipm_summary <- function(result){

  summ <- simanalyse::sma_summarise(result, measures=c("mean", "lower.q", "upper.q", "var"))

  p.mean <- summ$mcmcr1$mean$p
  f.mean <- summ$mcmcr1$mean$f
  N1.mean <- summ$mcmcr1$mean$N1
  sig.mean <- summ$mcmcr1$mean$sigma
  phi.mean <- summ$mcmcr1$mean$phi
  xi.mean <- summ$mcmcr1$mean$xi

  p.var <- summ$mcmcr1$var$p
  f.var <- summ$mcmcr1$var$f
  N1.var <- summ$mcmcr1$var$N1
  sig.var <- summ$mcmcr1$var$sigma
  phi.var <- summ$mcmcr1$var$phi
  xi.var <- summ$mcmcr1$var$xi

  p.lb <- summ$mcmcr1$lower.q$p
  f.lb <- summ$mcmcr1$lower.q$f
  N1.lb <- summ$mcmcr1$lower.q$N1
  sig.lb <- summ$mcmcr1$lower.q$sigma
  phi.lb <- summ$mcmcr1$lower.q$phi
  xi.lb <- summ$mcmcr1$lower.q$xi

  p.ub <- summ$mcmcr1$upper.q$p
  f.ub <- summ$mcmcr1$upper.q$f
  N1.ub <- summ$mcmcr1$upper.q$N1
  sig.ub <- summ$mcmcr1$upper.q$sigma
  phi.ub <- summ$mcmcr1$upper.q$phi
  xi.ub <- summ$mcmcr1$upper.q$xi

  res <- rbind(c(p.mean,f.mean,N1.mean,sig.mean,phi.mean,xi.mean),
               c(p.lb,f.lb,N1.lb,sig.lb,phi.lb,xi.lb),
               c(p.ub,f.ub,N1.ub,sig.ub,phi.ub,xi.ub),
               c(p.var,f.var,N1.var,sig.var,phi.var,xi.var))

  K <- length(f.mean)

  colnames(res) <- c(paste0("p[",1:(K-1),"]"), paste0("f[",1:K,"]"), "N1", "sigma",
                     paste0("phi[",1:(K-1),"]"), paste0("xi[",1:K,"]"))

  rownames(res) <- c("mean","var","2.5%","97.5%")

  return(res)
}



