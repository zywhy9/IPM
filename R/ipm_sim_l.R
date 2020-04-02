# Contents of file `ipm_sim_l.R`

#' Simulate population by the true joint likelihood model
#'
#' @param K a scalar, number of years.
#' @param N.1 a scalar,size of initial population.
#' @param sigma a scalar, standard deviation of count data.
#' @param p K x 1 vector, recapture probability of every year.
#' @param xi K x 1 vector, first-time capture probability of every year.
#' @param f K x 1 vector, fecundity rate of every year.
#' @param theta K x 2 matrix, mortality rate of every step.
#' @return data, the simulation resule
#' @export

ipm_sim_l <- function(K=4,
                      N.1=500,
                      sig=30,
                      p=rep(0.8,4),
                      xi=rep(0.8,4),
                      f=rep(2,4),
                      theta=matrix(rep(0.9,8),4,2)){

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
        C[K] ~ dbin(xi[K],Nu[K,2]+B[K])
        for(j in 1:(K-1)){
            M[j,K] ~ dbin(p[K],Sr[j,K,2])
        }
        Sr[K,K,2] <- C[K] + sum(M[1:(K-1),K])

        # Production
        phi[1] <- theta[1,1]*theta[1,2]
        phi[2] <- theta[1,2]*theta[2,1]
        for(i in 1:(K-2)){
            phi[2*i+1] <- theta[(i+1),1]*theta[(i+1),2]
            phi[2*(i+1)] <- theta[(i+1),2]*theta[(i+2),1]
        }
        phi[(2*K-1)] <- theta[K,1]*theta[K,2]
  "

  ## Simulation
  data <- sims_simulate(truemod, constants=nlist(K=K), latent=NA, parameters=nlist(N.1=N.1, sigma=sig, p=p, f=f, theta=theta, xi=xi))
  return(data)
}
