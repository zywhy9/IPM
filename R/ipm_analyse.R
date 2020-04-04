# Contents of file `ipm_analyse.R`

#' Analyse data using composite likelihood for a simulation study
#'
#' Using the composite likelihood model to analyse the simulation by true joint likelihood model.
#'
#' @param data an nlists object, simulation dataset returned by ipm_sim_l.
#' @param Plot a flag, indicates whether to save the traceplot and the density plot for MCMC outputs.
#' @param priors a string of code to set the prior.
#' @param maxtime a scalar, specifying the maximum time to spend on analysis.
#' @param unit a character string specifying the units of time for \code{max.time}. See \code{difftime}.
#' @param save a scalar, the number of (potentially thinned) samples to save.
#' @param chain a scalar, the number of MCMC chains.
#' @return a mcmcrs object, the MCMC outputs, and plots (If plot=TRUE).
#' @export


ipm_analyse <- function(data,
                        Plot=FALSE,
                        priors=NULL,
                        maxtime=5,
                        unit="mins",
                        save=20000L,
                        chain=3){

  if(is.null(priors)){
    priors <- "N.1 <- round(N1)
          N1 ~ dunif(0,2000)       # Initial total population
          sigma ~ dunif(0,100)     # Standard deviation of counting number
          for(i in 1:(K-1)){
              p[i] ~ dbeta(1,1)    # Recapture probability
              phi[i] ~ dbeta(1,1)  # Mortality rate
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
  compmod <- " N[1] <- 0
      for(i in 1:(K-1)){
          probas[i,i+1] <- phi[i]*p[i]
      }
      for(i in 1:(K-2)){
          for(j in (i+2):K){
             probas[i,j] <- prod(phi[i:(j-1)])*p[j-1]*prod(1-p[i:(j-2)])
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
      N[2] ~ dbin(phi[1],N.1+B[1])

      for(i in 3:K){
          Y[i] ~ dnorm(N[i],1/sigma/sigma)
          B[i-1] ~ dpois(N[i-1]*f[i-1]/2)
          N[i] ~ dbin(phi[i-1],N[i-1]+B[i-1])
      }

      B[K] ~ dpois(N[K]*f[K]/2)

      # Not Marked
      Nu[1] <- N.1
      for(i in 2:K){
          Nu[i] <- N[i] - sum(M[1:(i-1),i])/p[(i-1)]
      }

      # First Capture
      xi[1] <- R[1]/(Nu[1]+B[1])
      for(i in 2:K){
          xi[i] <- (R[i]-sum(M[1:(i-1),i]))/(Nu[i]+B[i])
      }

  "

  ## Transform
  K <- data[[1]]$K
  R <- rep(NA,K)
  Z <- rep(NA,(K-1))
  for(i in 1:(K-1)){
    Z[i] <- data[[1]]$Sr[i,i]-sum(data[[1]]$M[i,(i+1):K])
    R[i] <- data[[1]]$Sr[i,i]
  }
  R[K] <- data[[1]]$Sr[K,K]

  data[[1]]$M <- cbind(data[[1]]$M,Z)
  data <- subset(data, pars=c("M","Y","K"))
  data[[1]]$R <- R

  ## Result
  out <- simanalyse::sma_analyse_bayesian(data,
                                          compmod,
                                          priors,
                                          monitor=c("p","f","N1","sigma","phi","xi"),
                                          inits = list("N1"=500, f=rep(4,K),phi=rep(0.99,(K-1))),
                                          mode=sma_set_mode("paper", n.chains=chain, max.time=maxtime,
                                                            units=unit, n.save=save))
                                          # mode=sma_set_mode("quick"))
  out2 <- subset(mcmcr::as.mcmcrs(out), pars=c("p","f","N1","sigma","phi","xi"))

  if(Plot){
    pdf("plot.pdf")
    plot(window(coda::as.mcmc.list(out2[[1]])))
    graphics.off()
  }

  return(out2)
}



