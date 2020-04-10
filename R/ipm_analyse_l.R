# Contents of file `ipm_analyse_l.R`

#' Analyse data using composite likelihood for a simulation study
#'
#' Using the true joint likelihood model to analyse the simulation by true joint likelihood model.
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


ipm_analyse_l <- function(data,
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
              xi[i] ~ dbeta(1,1)   # First capture probability
              f[i] ~ dunif(0,10)   # Fecundity rate
          }
    "
  }

  ## True joint likelihood model
  truemod <- "
    # Year 1 Period 1 (Count Data)
    Nu[1] <- N.1                             # Nu[i] : Number of unmarked individuals in year i
    Nm[1] <- 0                               # Nm[i] : Number of marked individuals in year i
    Nt[1] <- Nu[1] + Nm[1]                   # Nt[i] : Number of total population in year i
    Y[1] ~ dnorm(Nt[1],1/(sigma^2))        # Y[i] : Count number of population in year i
    # Year 1 Period 2 (Birth & Capture-recapture)
    B[1] ~ dpois(Nt[1]*f[1]/2)               # B[i] : Number of newborn in year i
    C[1] ~ dbin(xi[1],Nu[1]+B[1])            # C[i] : Number of first captured individuals in year i
    Sr[1,1] <- C[1]                          # Sr[i,j]: Number of survived marked individuals who were captured in year i, in year j.

    # Year 2 Period 1
    Nu[2] ~ dbin(phi[1],Nu[1]+B[1]-C[1])
    Sr[1,2] ~ dbin(phi[1],Sr[1,1])
    Nm[2] <- Sr[1,2]
    Nt[2] <- Nu[2] + Nm[2]
    Y[2] ~ dnorm(Nt[2],1/(sigma^2))
    # Year 2 Period 2
    B[2] ~ dpois(Nt[2]*f[2]/2)
    C[2] ~ dbin(xi[2],Nu[2]+B[2])
    M[1,2] ~ dbin(p[1],Sr[1,2])              # M[i,j] : Number of individuals released in year i and recaptured in year j
    Sr[2,2] <- C[2] + M[1,2]
"
  model3 <- "# Year 3 Period 1
    Nu[3] ~ dbin(phi[2],Nu[2]+B[2]-C[2])
    Sr[1,3] ~ dbin(phi[2],Sr[1,2]-M[1,2])
    Sr[2,3] ~ dbin(phi[2],Sr[2,2])
    Nm[3] <- Sr[1,3] + Sr[2,3]
    Nt[3] <- Nu[3] + Nm[3]
    Y[3] ~ dnorm(Nt[3],1/(sigma^2))
    # Year 3 Period 2
    B[3] ~ dpois(Nt[3]*f[3]/2)
    C[3] ~ dbin(xi[3],Nu[3]+B[3])
    M[1,3] ~ dbin(p[2],Sr[1,3])
    M[2,3] ~ dbin(p[2],Sr[2,3])
    Sr[3,3] <- C[3] + sum(M[1:2,3])
"
  model4 <- "# Year K Period 1
      Nu[K] ~ dbin(phi[(K-1)],Nu[(K-1)]+B[(K-1)]-C[(K-1)])
      for(j in 1:(K-2)){
          Sr[j,K] ~ dbin(phi[(K-1)],Sr[j,(K-1)]-M[j,(K-1)])
      }
      Sr[(K-1),K] ~ dbin(phi[(K-1)],Sr[(K-1),(K-1)])
      Nm[K] <- sum(Sr[1:(K-1),K])
      Nt[K] <- Nu[K] + Nm[K]
      Y[K] ~ dnorm(Nt[K],1/(sigma^2))
      # Year K Period 2
      B[K] ~ dpois(Nt[K]*f[K]/2)
      C[K] ~ dbin(xi[K],Nu[K]+B[K])
      for(j in 1:(K-1)){
          M[j,K] ~ dbin(p[(K-1)],Sr[j,K])
      }
      Sr[K,K] <- C[K] + sum(M[1:(K-1),K])
"
  model4m <-  "# Year 4 to K-1
  for(i in 4:(K-1)){
      # Period 1
      Nu[i] ~ dbin(phi[(i-1)],Nu[(i-1)]+B[(i-1)]-C[(i-1)])
      for(j in 1:(i-2)){
          Sr[j,i] ~ dbin(phi[(i-1)],Sr[j,(i-1)]-M[j,(i-1)])
      }
      Sr[(i-1),i] ~ dbin(phi[(i-1)],Sr[(i-1),(i-1)])
      Nm[i] <- sum(Sr[1:(i-1),i])
      Nt[i] <- Nu[i] + Nm[i]
      Y[i] ~ dnorm(Nt[i],1/(sigma^2))
      # Period 2
      B[i] ~ dpois(Nt[i]*f[i]/2)
      C[i] ~ dbin(xi[i],Nu[i]+B[i])
      for(j in 1:(i-1)){
          M[j,i] ~ dbin(p[(i-1)],Sr[j,i])
      }
      Sr[i,i] <- C[i] + sum(M[1:(i-1),i])
  }
"

  ## Process
  K <- length(data[[1]]$B)
  if(K < 4){
    data[[1]]$K <- K
  }
  data <- subset(data, pars=c("C","M","Y","K"))

  ## Result
  if(K==3){
    truemod <- paste(truemod, model3, sep="\n")
  }else if(K==4){
    truemod <- paste(truemod, model3, model4, sep="\n")
  }else if(K>4){
    truemod <- paste(truemod, model3, model4, model4m, sep="\n")
  }

  out <- simanalyse::sma_analyse_bayesian(data,
                                          truemod,
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



