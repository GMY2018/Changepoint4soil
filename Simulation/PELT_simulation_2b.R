## --------------------------------
## PELT simulation scenario 2
## --------------------------------

source("PELT_spike_function_sim_nlm.R")


## Scenario 2: Seasonal patterns in the occurrence of the spikes
N <- 200
nt <- 5000
pen.sim <- 200
noise <- c(0.0005, 0.001, 0.005)

# Define the list to store the true values
spike.list <- vector(mode="list", length=N)
delta.list <- gamma.list <- asym.list <- vector(mode="list", length=N)
cpt.list <- vector(mode="list", length=N)
coef.list <- vector(mode="list", length=N)
mse.vec <- rep(0, times=N)
bad.vec <- timing.vec <- rep(0, times=N)

for (it in 1:N) {
  set.seed(it)
  
  np1 <- np2 <- 0
  while (np1 == 0) {
    spk.vec1 <- which(rpois(n=floor(nt/2), 0.002) == 1)
    spk.vec1 <- spk.vec1[spk.vec1 > 1]
    np1 <- length(spk.vec1)
  }
  while (np2 == 0) {
    spk.vec2 <- which(rpois(n=nt-floor(nt/2), 0.005) == 1) + floor(nt/2)
    spk.vec2 <- spk.vec2[spk.vec2 < nt]
    np2 <- length(spk.vec2)
  }
  spk.vec <- c(spk.vec1, spk.vec2 )
  
  delta.vec1 <- c(0, runif(n=np1, min=0.1, max=0.12))
  delta.vec2 <- runif(n=np2, min=0.05, max=0.1)
  delta.vec <- c(delta.vec1, delta.vec2)
  
  ns <- length(delta.vec)
  gamma.vec1 <- runif(n=length(delta.vec1), min=0.99, max=0.995)
  gamma.vec2 <- runif(n=length(delta.vec2), min=0.95, max=0.99)
  gamma.vec <- c(gamma.vec1, gamma.vec2)
  
  asym.vec <- runif(ns, min=0.05, max=0.08)
  
  sm.sim <- rep(NA, nt)
  sID <- c(1, spk.vec+1)
  eID <- c(spk.vec, nt)
  sgl <- eID[1] - sID[1]
  sm.sim[sID[1]:eID[1]] <- (runif(n=1, 0.1, 0.2) + delta.vec[1] - asym.vec[1])*(gamma.vec[1]^(0:sgl)) + 
    asym.vec[1] + rnorm(sgl+1, mean=0, sd=noise[2])
  for (i in 2:length(sID)) {
    sgl <- eID[i] - sID[i]
    sm.last <- sm.sim[eID[i-1]]
    sm.sim[sID[i]:eID[i]] <- (sm.last + delta.vec[i] - asym.vec[i])*(gamma.vec[i]^(0:sgl)) + 
      asym.vec[i] + rnorm(sgl+1, mean=0, sd=noise[2])
  }
  
  # Apply the algorithm
  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl2(data=sm.sim, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.sim, minsl=24, costfun=spike.exp.nlm1, 
                            upper.par=c(0.5, 0.7, 1), thresh=0.003, nprune=FALSE)
  sm.fit <- spike.PELT.yhat(cpt=c(0, sm.pelt$cpt), data=sm.sim, xreg=0,
                            lastchangecoef=sm.pelt$lastchangecoef, type=1)
  time2 <- proc.time() - time1
  
  # Record the result
  spike.list[[it]] <- spk.vec
  delta.list[[it]] <- delta.vec
  gamma.list[[it]] <- gamma.vec
  asym.list[[it]] <- asym.vec 
  cpt.list[[it]] <- sm.pelt$cpt 
  coef.list[[it]] <- sm.pelt$lastchangecoef[sm.pelt$cpt+1]  
  
  mse.vec[it] <- mean((sm.sim - sm.fit$sm.hat)^2, na.rm=TRUE)
  bad.vec[it] <- sum(is.na(sm.fit$sm.hat))
  timing.vec[it] <- time2[1]
  
  if (it %% 10 == 0) {
    save(spike.list, delta.list, gamma.list, asym.list, cpt.list, coef.list,
         mse.vec, bad.vec, timing.vec, file="PELT-sim-nlm-s2b.RData")
  }
}

save(spike.list, delta.list, gamma.list, asym.list, cpt.list, coef.list,
     mse.vec, bad.vec, timing.vec, file="PELT-sim-nlm-s2b.RData")

