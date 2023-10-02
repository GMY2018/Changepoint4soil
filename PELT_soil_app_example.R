## -------------------------------------------
## NEON soil moisture time series analysis
## -------------------------------------------

## Changepoint analysis of soil moisture time series from
## field site SRER, SERC and TALL

source("PELT_spike_function_nlm.R")


## Field site SERC
load("PELT_SERC_data_loc1.RData")
sm.dat2 <- sm.dat1[seq(2, length(time0), by=2)]   # hourly data
nt <- length(sm.dat2)

## Run the model using different penalty parameters
pen.list <- c(100, 150, 200, 250, 300)
np <- length(pen.list)

## Fitting using the new function
PELT.result.loc1 <- vector(mode="list", length=np)
timing.loc1 <- rep(0, np)

# Model: y_t = asym + alpha0*exp(-exp(lgamma)*t) + e_t
for (i in 1:np) {
  # choose a penalty
  pen.i <- pen.list[i]

  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=sm.dat2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=12, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 3), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  PELT.result.loc1[[i]] <- sm.pelt
  timing.loc1[i] <- time2[1]

  print(i)

  save(PELT.result.loc1, timing.loc1, file="PELT_SERC_app_1.RData")
}


## Field site SRER
load("PELT_SRER_data_loc4.RData")
sm.dat2 <- sm.dat1[seq(2, length(time0), by=2)]
nt <- length(sm.dat2)

## Run the model using different penalty parameters
pen.list <- c(50, 100, 150, 200, 250)
np <- length(pen.list)

## Fitting using the new function
PELT.result.loc4 <- vector(mode="list", length=np)
timing.loc4 <- rep(0, np)

# Model: y_t = asym + alpha0*exp(-exp(lgamma)*t) + e_t
for (i in 1:np) {
  # choose a penalty
  pen.i <- pen.list[i]

  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=sm.dat2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=12, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 3), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  PELT.result.loc4[[i]] <- sm.pelt
  timing.loc4[i] <- time2[1]

  print(i)

  save(PELT.result.loc4, timing.loc4, file="PELT_SRER_app_1.RData")
}


## Field site TALL
load("PELT_TALL_data_loc5.RData")
tID <- (time0 >= as.POSIXct("2018-02-01", tz="UTC")) & (time0 <= as.POSIXct("2019-02-01", tz="UTC"))
time0 <- time0[tID]
sm.dat1 <- sm.dat1[tID]
sm.dat2 <- sm.dat1[seq(2, length(time0), by=2)]
nt <- length(sm.dat2)

## Run the model using different penalty parameters
pen.list <- c(50, 100, 150, 200, 250)
np <- length(pen.list)

## Fitting using the new function
PELT.result.loc5 <- vector(mode="list", length=np)
timing.loc5 <- rep(0, np)

# Model: y_t = asym + alpha0*exp(-exp(lgamma)*t) + e_t
for (i in 1:np) {
  # choose a penalty
  pen.i <- pen.list[i]

  time1 <- proc.time()
  sm.pelt <- spike.PELT.msl(data=sm.dat2, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=12, costfun=spike.exp.nlm2,
                            upper.par=c(0.4, 0.5, 3), thresh=0.001, nprune=FALSE)
  time2 <- proc.time() - time1
  PELT.result.loc5[[i]] <- sm.pelt
  timing.loc5[i] <- time2[1]

  print(i)

  save(PELT.result.loc5, timing.loc5, file="PELT_TALL_app_1.RData")
}


