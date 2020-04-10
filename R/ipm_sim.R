# Contents of file `ipm_sim.R`

#' Simulation
#'
#' Simulate population by the model of interest.
#'
#' @param K a scalar, number of years.
#' @param N.1 a scalar,size of initial population.
#' @param sigma a scalar, standard deviation of count data.
#' @param p K x 1 vector, recapture probability of every year.
#' @param xi K x 1 vector, first-time capture probability of every year.
#' @param f K x 1 vector, fecundity rate of every year.
#' @param phi (K-1) x 1 vector, mortality rate of every year.
#' @param R (K-1) x 1 vector, number of release in every year.
#' @param model a string, specifying the model of interest.
#' @param param an nlist object, specifying the parameters of interest.
#' @return data, the simulation result.
#' @export

ipm_sim <- function(K=4,
                    N.1=500,
                    sigma=30,
                    p=rep(0.8,4),
                    xi=rep(0.8,4),
                    f=rep(2,4),
                    phi=rep(0.9,3),
                    R=NULL,
                    model="l",
                    param=NULL){
  if(model=="l"){
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

    if(K==3){
      truemod <- paste(truemod, model3, sep="\n")
    }else if(K==4){
      truemod <- paste(truemod, model3, model4, sep="\n")
    }else if(K>4){
      truemod <- paste(truemod, model3, model4, model4m, sep="\n")
    }

    if(length(p)==K){
      p <- p[1:(K-1)]
    }

    ## Simulation
    if(K>3){
      data <- sims::sims_simulate(truemod, constants=nlist::nlist(K=K), latent=NA,
                                  parameters=nlist::nlist(N.1=N.1, sigma=sigma, p=p, f=f, phi=phi, xi=xi))
    }else{
      data <- sims::sims_simulate(truemod, latent=NA,
                                  parameters=nlist::nlist(N.1=N.1, sigma=sigma, p=p, f=f, phi=phi, xi=xi))
    }
  }else if(model=="cl"){
    if(is.null(R)){
      R <- c(N.1*1.5, N.1*2.5, N.1*4)
    }

    ## Composite likelihood model
    compmod <- " N[1] <- 0
      for(i in 1:(K-1)){
          probas[i,i+1] <- phi[i]*p[i]
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

      B[K] ~ dpois(N[K]*f[K]/2)
    "
    model3m <- "
      for(i in 1:(K-2)){
          for(j in (i+2):K){
             probas[i,j] <- prod(phi[i:(j-1)])*p[j-1]*prod(1-p[i:(j-2)])
          }
      }

       for(i in 3:K){
          Y[i] ~ dnorm(N[i],1/sigma/sigma)
          B[i-1] ~ dpois(N[i-1]*f[i-1]/2)
          N[i] ~ dbin(phi[i-1],N[i-1]+B[i-1])
      }
    "
    if(K>2){
      compmod <- paste(compmod, model3m, sep="\n")
    }



    ## Simulation
    data <- sims::sims_simulate(compmod, constants=nlist::nlist(K=K), latent=NA,
                                parameters=nlist::nlist(N.1=N.1, sigma=sigma, p=p, f=f, phi=phi, R=R))
  }else{
    data <- sims::sims_simulate(model,latent=NA,parameters=param)
  }

  return(data)
}
