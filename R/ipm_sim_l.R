# Contents of file `ipm_sim_l.R`

#' True joint likelihood simulation
#'
#' Simulate population by the true joint likelihood model
#'
#' @param K a scalar, number of years.
#' @param N.1 a scalar,size of initial population.
#' @param sigma a scalar, standard deviation of count data.
#' @param p K x 1 vector, recapture probability of every year.
#' @param xi K x 1 vector, first-time capture probability of every year.
#' @param f K x 1 vector, fecundity rate of every year.
#' @param phi (K-1) x 1 vector, mortality rate of every year.
#' @return data, the simulation resule
#' @export

ipm_sim_l <- function(K=4,
                      N.1=500,
                      sigma=30,
                      p=rep(0.8,3),
                      xi=rep(0.8,4),
                      f=rep(2,4),
                      phi=rep(0.9,3)){

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

        # Year 3 Period 1
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

        # Year 4 to K-1
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

        # Year K Period 1
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

  ## Simulation
  data <- sims::sims_simulate(truemod, constants=nlist::nlist(K=K), latent=NA,
                              parameters=nlist::nlist(N.1=N.1, sigma=sigma, p=p, f=f, phi=phi, xi=xi))
  return(data)
}
