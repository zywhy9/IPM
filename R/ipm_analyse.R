# Contents of file `ipm_analyse.R`

#' Using the true joint likelihood model to analyse the simulation by composite likelihood model
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
          N1 ~ dunif(0,2000)  # Initial total population
          sigma ~ dunif(0,100)       # Standard deviation of counting number
          for(i in 1:K){
              p[i] ~ dbeta(1,1)    # Recapture probability
              xi[i] ~ dbeta(1,1)   # First capture probability
              f[i] ~ dunif(0,10)   # Fecundity rate
              for(j in 1:2){
                theta[i,j] ~ dunif(0,1)
              }
          }
    "
  }

  ## True joint likelihood model
  truemod <- "
        # Year 1 Period 1 (Count Data)
        Nu[1,1] <- N.1                           # Nu[i,j] : Number of unmarked individuals in year i period j
        Nm[1,1] <- 0                             # Nm[i,j] : Number of marked individuals in year i period j
        Nt[1,1] <- Nu[1,1] + Nm[1,1]             # Nt[i,j] : Number of total population in year i period j
        Y[1] ~ dnorm(Nt[1,1],1/(sigma^2))        # Y[i] : Count number of population in year i
        # Year 1 Period 2 (Birth & Capture-recapture)
        Nu[1,2] ~ dbin(theta[1,1],Nu[1,1])
        Nt[1,2] <- Nu[1,2]
        B[1] ~ dpois(Nt[1,2]*f[1]/2)             # B[i] : Number of newborn in year i
        C[1] ~ dbin(xi[1],Nu[1,2]+B[1])          # C[i] : Number of first captured individuals in year i
        Sr[1,1,2] <- C[1]                        # Sr[i,j,k]: Number of survived marked individuals who were captured in year i, at time point year j period k.
        Nm[1,2] <- Sr[1,1,2]

        # Year 2 Period 1
        Nu[2,1] ~ dbin(theta[1,2],Nu[1,2]+B[1]-C[1])
        Sr[1,2,1] ~ dbin(theta[1,2],Sr[1,1,2])
        Nm[2,1] <- Sr[1,2,1]
        Nt[2,1] <- Nu[2,1] + Nm[2,1]
        Y[2] ~ dnorm(Nt[2,1],1/(sigma^2))
        # Year 2 Period 2
        Nu[2,2] ~ dbin(theta[2,1],Nu[2,1])
        Sr[1,2,2] ~ dbin(theta[2,1],Sr[1,2,1])
        Nm[2,2] <- Sr[1,2,2]
        Nt[2,2] <- Nu[2,2] + Nm[2,2]
        B[2] ~ dpois(Nt[2,2]*f[2]/2)
        C[2] ~ dbin(xi[2],Nu[2,2]+B[2])
        M[1,2] ~ dbin(p[2],Sr[1,2,2])            # M[i,j] : Number of individuals released in year i and recaptured in year j
        Sr[2,2,2] <- C[2] + M[1,2]

        # Year 3 Period 1
        Nu[3,1] ~ dbin(theta[2,2],Nu[2,2]+B[2]-C[2])
        Sr[1,3,1] ~ dbin(theta[2,2],Sr[1,2,2]-M[1,2])
        Sr[2,3,1] ~ dbin(theta[2,2],Sr[2,2,2])
        Nm[3,1] <- Sr[1,3,1] + Sr[2,3,1]
        Nt[3,1] <- Nu[3,1] + Nm[3,1]
        Y[3] ~ dnorm(Nt[3,1],1/(sigma^2))
        # Year 3 Period 2
        Nu[3,2] ~ dbin(theta[3,1],Nu[3,1])
        Sr[1,3,2] ~ dbin(theta[3,1],Sr[1,3,1])
        Sr[2,3,2] ~ dbin(theta[3,1],Sr[2,3,1])
        Nm[3,2] <- Sr[1,3,2] + Sr[2,3,2]
        Nt[3,2] <- Nu[3,1] + Nm[3,1]
        B[3] ~ dpois(Nt[3,2]*f[3]/2)
        C[3] ~ dbin(xi[3],Nu[3,2]+B[3])
        M[1,3] ~ dbin(p[3],Sr[1,3,2])
        M[2,3] ~ dbin(p[3],Sr[2,3,2])
        Sr[3,3,2] <- C[3] + sum(M[1:2,3])

        # Year 4 to K-1
        for(i in 4:(K-1)){
            # Period 1
            Nu[i,1] ~ dbin(theta[(i-1),2],Nu[(i-1),2]+B[(i-1)]-C[(i-1)])
            for(j in 1:(i-2)){
                Sr[j,i,1] ~ dbin(theta[(i-1),2],Sr[j,(i-1),2]-M[j,(i-1)])
            }
            Sr[(i-1),i,1] ~ dbin(theta[(i-1),2],Sr[(i-1),(i-1),2])
            Nm[i,1] <- sum(Sr[1:(i-1),i,1])
            Nt[i,1] <- Nu[i,1] + Nm[i,1]
            Y[i] ~ dnorm(Nt[i,1],1/(sigma^2))
            # Period 2
            Nu[i,2] ~ dbin(theta[i,1],Nu[i,1])
            for(j in 1:(i-1)){
                Sr[j,i,2] ~ dbin(theta[i,1],Sr[j,i,1])
            }
            Nm[i,2] <- sum(Sr[1:(i-1),i,2])
            Nt[i,2] <- Nu[i,2] + Nm[i,2]
            B[i] ~ dpois(Nt[i,2]*f[i]/2)
            C[i] ~ dbin(xi[i],Nu[i,2]+B[i])
            for(j in 1:(i-1)){
                M[j,i] ~ dbin(p[i],Sr[j,i,2])
            }
            Sr[i,i,2] <- C[i] + sum(M[1:(i-1),i])
        }

        # Year K Period 1
        Nu[K,1] ~ dbin(theta[(K-1),2],Nu[(K-1),2]+B[(K-1)]-C[(K-1)])
        for(j in 1:(K-2)){
            Sr[j,K,1] ~ dbin(theta[(K-1),2],Sr[j,(K-1),2]-M[j,(K-1)])
        }
        Sr[(K-1),K,1] ~ dbin(theta[(K-1),2],Sr[(K-1),(K-1),2])
        Nm[K,1] <- sum(Sr[1:(K-1),K,1])
        Nt[K,1] <- Nu[K,1] + Nm[K,1]
        Y[K] ~ dnorm(Nt[K,1],1/(sigma^2))
        # Year K Period 2
        Nu[K,2] <- Nu[K,1] - (Nu[(K-1),2] + B[(K-1)] - C[(K-1)] - Nu[K,1])
        for(j in 1:(K-2)){
            Sr[j,K,2] <- Sr[j,K,1] - (Sr[j,(K-1),2] - M[j,(K-1)] - Sr[j,K,1])
        }
        Sr[(K-1),K,2] <- Sr[(K-1),K,1] - (Sr[(K-1),(K-1),2] - Sr[(K-1),K,1])
        Nm[K,2] <- sum(Sr[1:(K-1),K,2])
        Nt[K,2] <- Nu[K,2] + Nm[K,2]
        B[K] ~ dpois(Nt[K,2]*f[K]/2)
        for(j in 1:(K-1)){
            M[j,K] ~ dbin(p[K],Sr[j,K,2])
        }

        # Production
        phi[1] <- theta[1,1]*theta[1,2]
        phi[2] <- theta[1,2]*theta[2,1]
        for(i in 1:(K-2)){
            phi[2*i+1] <- theta[(i+1),1]*theta[(i+1),2]
            phi[2*(i+1)] <- theta[(i+1),2]*theta[(i+2),1]
        }
        phi[(2*K-1)] <- theta[K,1]*theta[K,2]
  "

  ## Transform
  K <- data[[1]]$K
  M <- ifelse(is.na(data[[1]]$M),0,data[[1]]$M)
  C <- apply(M,1,sum)
  M <- M[,1:K]
  data <- nlists(nlist(M=M,K=K,Y=data[[1]]$Y,C=C))

  ## Result
  out <- simanalyse::sma_analyse_bayesian(data,
                                          truemod,
                                          priors,
                                          monitor=c("p","f","N1","sigma","theta","phi"),
                                          inits = list("N1"=500, f=rep(4,K),theta=matrix(rep(0.99,K*2),K,2)),
                                          mode=sma_set_mode("paper", n.chains=chain, max.time=maxtime, units=unit, n.save=save))
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
  colnames(res) <- c("p[1]","p[2]","p[3]","p[4]","f[1]","f[2]","f[3]","f[4]","N1","sigma",
                     "phi[1]","phi[2]","phi[3]","phi[4]","phi[5]","phi[6]","phi[7]","theta[1,1]",
                     "theta[1,2]","theta[2,1]","theta[2,2]","theta[3,1]","theta[3,2]","theta[4,1]",
                     "theta[4,2]")
  rownames(res) <- c("mean","2.5%","97.5%")

  return(res)
}



