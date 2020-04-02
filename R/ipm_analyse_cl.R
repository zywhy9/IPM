# Contents of file `ipm_analyse_cl.R`

#' Analysing the simulation by composite likelihood model
#'
#' @param data an nlist object, simulation dataset returned by ipm_sim_cl.
#' @param Plot a flag, indicates whether to save the traceplot and the density plot for MCMC outputs.
#' @param prior a string of code to set the prior.
#' @param maxtime a scalar, specifying the maximum time to spend on analysis.
#' @param unit a character string specifying the units of time for \code{max.time}. See \code{difftime}.
#' @param save a scalar, the number of (potentially thinned) samples to save.
#' @param chain a scalar, the number of MCMC chains.
#' @return a matirx, with summary staitstics, and plots (If plot=TRUE).
#' @export


ipm_analyse <- function(data,
                        Plot=FALSE,
                        prior=NULL,
                        maxtime=10,
                        unit="mins",
                        save=20000L,
                        chain=3){

  if(is.null(prior)){
    priors <- "N.1 <- round(N1)
          N1 ~ dunif(0,2000)       # Initial total population
          sigma ~ dunif(0,100)     # Standard deviation of counting number
          for(i in 1:(K-1)){
              p[i] ~ dbeta(1,1)    # Recapture probability
          }
          for(i in 1:K){
              f[i] ~ dunif(0,10)   # Fecundity rate
              for(j in 1:2){
                  theta[i,j] ~ dunif(0,1)
              }
          }
    "
  }

  ## Composite likelihood model
  compmod <- "
      for(i in 1:(K-1)){
          probas[i,i+1] <- theta[i,2]*p[i]*theta[(i+1),1]
      }
      for(i in 1:(K-2)){
          for(j in (i+2):K){
             probas[i,j] <- prod(theta[i:(j-1),2]*theta[(i+1):j,1])*p[j-1]*prod(1-p[i:(j-2)])
          }
      }

      for(i in 1:(K-1)){
          probas[i,K+1] <- 1-sum(probas[i,(i+1):K])
      }

      for(i in 1:(K-1)){
         M[i,(i+1):(K+1)] ~ dmulti(probas[i,(i+1):(K+1)],R[i])
      }

      Y[1] ~ dnorm(N.1,1/sigma/sigma)
      Y[2] ~ dnorm(N[2],1/sigma/sigma)
      B[1] ~ dpois(N.1*f[1]/2)
      N[2] ~ dbin(prod(theta[1,1:2]),N.1+B[1])

      for(i in 3:K){
          Y[i] ~ dnorm(N[i],1/sigma/sigma)
          B[i-1] ~ dpois(N[i-1]*f[i-1]/2)
          N[i] ~ dbin(prod(theta[(i-1),1:2]),N[i-1]+B[i-1])
      }

      B[K] ~ dpois(N[K]*f[K]/2)

      # Production
      for(i in 1:K){
          phi[i] <- prod(theta[i,1:2])
      }
  "

  ## Process
  M <- ifelse(is.na(data[[1]]$M),0,data[[1]]$M)
  R <- apply(M,1,sum)
  data <- nlists(nlist(B=data[[1]]$B,M=data[[1]]$M,K=data[[1]]$K,Y=data[[1]]$Y,R=R,N=data[[1]]$N))

  ## Result
  out <- simanalyse::sma_analyse_bayesian(data,
                                          compmod,
                                          priors,
                                          monitor=c("p","f","N1","sigma","theta","phi"),
                                          inits = list("N1"=500, f=rep(4,K),theta=matrix(rep(0.99,K*2),K,2)),
                                          mode=sma_set_mode("paper",  max.time=maxtime, units=unit, n.save=save))
  out2 <- subset(mcmcr::as.mcmcrs(out), pars=c("p","f","N1","sigma","theta","phi"))

  if(Plot){
    pdf("plot.pdf")
    plot(window(as.mcmc.list(out2[[1]])))
    graphics.off()
  }

  ## Output
  summ <- sma_summarise(out2)

  p.mean <- summ$mcmcr1$mean$p
  f.mean <- summ$mcmcr1$mean$f
  N1.mean <- summ$mcmcr1$mean$N1
  sig.mean <- summ$mcmcr1$mean$sigma
  phi.mean <- summ$mcmcr1$mean$phi
  theta.mean <- as.vector(t(summ$mcmcr1$mean$theta))

  p.lb <- summ$mcmcr1$lower.q$p
  f.lb <- summ$mcmcr1$lower.q$f
  N1.lb <- summ$mcmcr1$lower.q$N1
  sig.lb <- summ$mcmcr1$lower.q$sigma
  phi.lb <- summ$mcmcr1$lower.q$phi
  theta.lb <- as.vector(t(summ$mcmcr1$lower.q$theta))

  p.ub <- summ$mcmcr1$upper.q$p
  f.ub <- summ$mcmcr1$upper.q$f
  N1.ub <- summ$mcmcr1$upper.q$N1
  sig.ub <- summ$mcmcr1$upper.q$sigma
  phi.ub <- summ$mcmcr1$upper.q$phi
  theta.ub <- as.vector(t(summ$mcmcr1$upper.q$theta))

  res <- rbind(c(p.mean,f.mean,N1.mean,sig.mean,phi.mean,theta.mean),
               c(p.lb,f.lb,N1.lb,sig.lb,phi.lb,theta.lb),
               c(p.ub,f.ub,N1.ub,sig.ub,phi.ub,theta.ub))
  colnames(res) <- c("p[1]","p[2]","p[3]","f[1]","f[2]","f[3]","f[4]","N1","sigma","phi[1]",
                     "phi[2]","phi[3]","phi[4]","theta[1,1]", "theta[1,2]","theta[2,1]",
                    "theta[2,2]","theta[3,1]","theta[3,2]","theta[4,1]","theta[4,2]")
  rownames(res) <- c("mean","2.5%","97.5%")

  return(res)
}



