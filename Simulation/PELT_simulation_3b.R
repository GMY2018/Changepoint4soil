## --------------------------------
## PELT simulation scenario 3
## --------------------------------

source("PELT_spike_function_sim_nlm.R")

## Uncomment corresponding lines to run the algorithm with large or small penalty

## Scenario 3: Combining different resolutions (heavy + small rainfall)
N <- 100
nt <- 5000
pen.sim1 <- 100   
pen.sim2 <- 800
# 100 for small scale, 800 for large scale
noise <- c(0.0005, 0.001)


# Define the list to store the true values
spike.list <- vector(mode="list", length=N)
delta.list <- gamma.list <- asym.list <- vector(mode="list", length=N)
cpt1.list <- cpt2.list <- vector(mode="list", length=N)
coef1.list <- coef2.list <- vector(mode="list", length=N)
mse1.vec <- mse2.vec <- rep(0, times=N)
bad1.vec <- bad2.vec <- timing1.vec <- timing2.vec <- rep(0, times=N)

it <- 1
ct <- 1
while (it <= N) {
  set.seed(ct+10)
  
  # leave a long drying period
  nt2 <- 48*15
  t2 <- sample(1000:3500, size=1)
  spk.vec1 <- which(rpois(n=t2, 0.002) == 1)
  spk.vec2 <- which(rpois(n=nt-nt2-t2, 0.002) == 1) + t2 + nt2
  spk.vec1 <- spk.vec1[spk.vec1 > 1]
  spk.vec2 <- spk.vec2[spk.vec2 < nt]
  spk.vec <- c(spk.vec1, spk.vec2)
  np1 <- length(spk.vec1)
  np2 <- length(spk.vec2)
  
  if (np1 + np2 > 0) {
    delta.vec <- c(0, runif(n=np1+np2, min=0.1, max=0.12))
    
    if (np1 > 0) {
      gamma.vec1 <- runif(n=np1, min=0.98, max=0.99)
    } else {
      gamma.vec1 <- NULL
    }
    if (np2 > 0) {
      gamma.vec2 <- runif(n=np2, min=0.98, max=0.99)
    } else {
      gamma.vec2 <- NULL
    }
    gamma.vec <- c(gamma.vec1, 0.995, gamma.vec2)
    
    asym.vec <- runif(length(gamma.vec), min=0.05, max=0.08)
    
    sm.sim.a <- rep(NA, nt)
    sID <- c(1, spk.vec+1)
    eID <- c(spk.vec, nt)
    sgl <- eID[1] - sID[1]
    sm.sim.a[sID[1]:eID[1]] <- (runif(n=1, 0.1, 0.2) + delta.vec[1] - asym.vec[1])*(gamma.vec[1]^(0:sgl)) + 
      asym.vec[1] + rnorm(sgl+1, mean=0, sd=noise[2])
    for (i in 2:length(sID)) {
      sgl <- eID[i] - sID[i]
      sm.last <- sm.sim.a[eID[i-1]]
      sm.sim.a[sID[i]:eID[i]] <- (sm.last + delta.vec[i] - asym.vec[i])*(gamma.vec[i]^(0:sgl)) + 
        asym.vec[i] + rnorm(sgl+1, mean=0, sd=noise[2])
    }
    
    # Adding the smaller ups and downs
    if (np1 == 0) spk.vec1 <- 0
    if (np2 == 0) spk.vec2 <- nt
    nt2 <- spk.vec2[1] - rev(spk.vec1)[1]
    # again, make sure there is at least one changepoint
    np3 <- 0
    while (np3 == 0) {
      spk.vec3 <- which(rpois(n=nt2, 0.01) == 1)
      spk.vec3 <- spk.vec3[(spk.vec3 > 1) & (spk.vec3 < nt2)]
      np3 <- length(spk.vec3)
    }
    
    delta.vec3 <- c(0, runif(n=length(spk.vec3), min=0.01, max=0.02))
    gamma.vec3 <- runif(n=length(delta.vec3), min=0.95, max=0.99)
    
    sm.sim3 <- rep(NA, nt2)
    sID <- c(1, spk.vec3+1)
    eID <- c(spk.vec3, nt2)
    sgl <- eID[1] - sID[1]
    sm.sim3[sID[1]:eID[1]] <- (0.05 + delta.vec3[1])*(gamma.vec3[1]^(0:sgl))
    for (i in 2:length(sID)) {
      sgl <- eID[i] - sID[i]
      sm.last <- sm.sim3[eID[i-1]]
      sm.sim3[sID[i]:eID[i]] <- (sm.last + delta.vec3[i])*(gamma.vec3[i]^(0:sgl))
    }
    
    # combine the two series
    sm.sim.b <- rep(0, times=length(sm.sim.a))
    if (rev(spk.vec1)[1] > 1) {
      sm.sim.b[(rev(spk.vec1)[1]+1):spk.vec2[1]] <- sm.sim3
    } else {
      sm.sim.b[1:spk.vec2[1]] <- sm.sim3
    }
    sm.sim <- sm.sim.a + sm.sim.b
    
    # Apply the algorithm
    # small penalty
    # time1 <- proc.time()
    # sm.pelt1 <- spike.PELT.msl2(data=sm.sim, xreg=0, ini.par=NULL, ini.asym=0.01,
    #                             pen=pen.sim1, minsl=24, costfun=spike.exp.nlm1, 
    #                             upper.par=c(0.5, 0.7, 1), thresh=0.003, nprune=FALSE)
    # sm.fit1 <- spike.PELT.yhat(cpt=c(0, sm.pelt1$cpt), data=sm.sim, xreg=0,
    #                            lastchangecoef=sm.pelt1$lastchangecoef, type=1)
    # time2 <- proc.time() - time1
    
    # large penalty
    time3 <- proc.time()
    sm.pelt2 <- spike.PELT.msl2(data=sm.sim, xreg=0, ini.par=NULL, ini.asym=0.01,
                                pen=pen.sim2, minsl=24, costfun=spike.exp.nlm1,
                                upper.par=c(0.5, 0.7, 1), thresh=0.003, nprune=FALSE)
    sm.fit2 <- spike.PELT.yhat(cpt=c(0, sm.pelt2$cpt), data=sm.sim, xreg=0,
                               lastchangecoef=sm.pelt2$lastchangecoef, type=1)
    time4 <- proc.time() - time3
    
    # Record the result
    spike.list[[it]] <- list(large=spk.vec, small=spk.vec3 + rev(spk.vec1)[1])
    delta.list[[it]] <- list(large=delta.vec, small=delta.vec3)
    gamma.list[[it]] <- list(large=gamma.vec, small=gamma.vec3)
    asym.list[[it]] <- asym.vec 
    # cpt1.list[[it]] <- sm.pelt1$cpt 
    # coef1.list[[it]] <- sm.pelt1$lastchangecoef[sm.pelt1$cpt+1]  
    cpt2.list[[it]] <- sm.pelt2$cpt
    coef2.list[[it]] <- sm.pelt2$lastchangecoef[sm.pelt2$cpt+1]
    
    # mse1.vec[it] <- mean((sm.sim - sm.fit1$sm.hat)^2, na.rm=TRUE)
    # bad1.vec[it] <- sum(is.na(sm.fit1$sm.hat))
    # timing1.vec[it] <- time2[1]
    mse2.vec[it] <- mean((sm.sim - sm.fit2$sm.hat)^2, na.rm=TRUE)
    bad2.vec[it] <- sum(is.na(sm.fit2$sm.hat))
    timing2.vec[it] <- time4[1]
    
    it <- it + 1
  } else {
    # if np1 + np2 = 0, then do not count
    it <- it + 0
  }
  
  # but always increase the count
  ct <- ct + 1 
  
  if (it %% 5 == 0) {
    save(spike.list, delta.list, gamma.list, asym.list, 
         # cpt1.list, coef1.list, 
         cpt2.list, coef2.list,
         # mse1.vec, bad1.vec, timing1.vec,
         mse2.vec, bad2.vec, timing2.vec,
         file="PELT-sim-nlm-s3b-p2.RData")
  }
}
  
save(spike.list, delta.list, gamma.list, asym.list, 
     # cpt1.list, coef1.list, 
     cpt2.list, coef2.list,
     # mse1.vec, bad1.vec, timing1.vec,
     mse2.vec, bad2.vec, timing2.vec,
     file="PELT-sim-nlm-s3b-p2.RData")

