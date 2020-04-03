# Contents of file `ipm_sim_cl.R`

#' Simulate population by the composite likelihood model
#'
#' @param K a scalar, number of years.
#' @param N.1 a scalar,size of initial population.
#' @param sigma a scalar, standard deviation of count data.
#' @param p K x 1 vector, recapture probability of every year.
#' @param f K x 1 vector, fecundity rate of every year.
#' @param phi (K-1) x 1 vector, mortality rate of every year.
#' @param R (K-1) x 1 vector, number of release in every year.
#' @return data, the simulation result.
#' @export

ipm_sim_cl <- function(K=4,
                       N.1=500,
                       sigma=30,
                       p=rep(0.8,4),
                       f=rep(2,4),
                       phi=rep(0.9,3),
                       R=NULL){

  if(is.null(R)){
    R <- c(N.1*1.5, N.1*2.5, N.1*4)
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
  "

  ## Simulation
  data <- sims_simulate(compmod, constants=nlist(K=K), latent=NA, parameters=nlist(N.1=N.1, sigma=sigma, p=p, f=f, phi=phi, R=R))
  return(data)
}
