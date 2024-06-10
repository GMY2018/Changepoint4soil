## -----------------------------------------
## Penalty learning using rainfall data
## -----------------------------------------


source("PELT_soil_function_nlm.R")
load("DataExample/Rain_SM_data_SRER.RData")


## Prepare the data (pick a short segent to experiment on)
tID1 <- neon.rain$time %in% neon.sm.top$time
tID2 <- neon.sm.top$time %in% neon.rain$time
neon.rain <- neon.rain[tID1, ]
neon.sm.top <- neon.sm.top[tID2, ]

nts <- nrow(neon.sm.top)
if (nts %% 2 == 1) {
  neon.sm.top <- neon.sm.top[-nts, ]
  neon.rain <- neon.rain[-nts, ]
}
sm0 <- neon.sm.top$hloc_2
sm1 <- sm0[seq(1, length(sm0), by=2)]
rain0 <- neon.rain$precip
rain1 <- rain0[seq(1, length(rain0), by=2)] + rain0[seq(2, length(rain0), by=2)]
time1 <- neon.sm.top$time[seq(1, length(sm0), by=2)]

tID <- 8500:10000
plot(time1[tID], rain1[tID], type="h", xlab="Time", ylab="Precipitation")
plot(time1[tID], sm1[tID], type="l", xlab="Time", ylab="Soil water content")

time1 <- time1[tID]
sm1 <- sm1[tID]
rain1 <- rain1[tID]

# fill the small missing gaps in soil moisture
library(zoo)
ar1 <- ar(sm1, order.max=1, na.action=na.pass)$ar
rle.mis <- rle(is.na(sm1))
rle.tab <- data.frame(values=rle.mis$values, length=rle.mis$lengths,
                      end=cumsum(rle.mis$lengths),
                      start=c(1, cumsum(rle.mis$lengths)[-length(rle.mis$lengths)]+1))
jit.tab <- rle.tab[(rle.tab$values == 1) & (rle.tab$length > 10), ]
sm.dat0 <- sm.dat <- na.approx(sm1)
for (i in 1:nrow(jit.tab)) {
  sm.dat[jit.tab$start[i]:jit.tab$end[i]] <- jitter(sm.dat0[jit.tab$start[i]:jit.tab$end[i]], factor=3)
}


## A sequence of penalty parameters
pen.list <- seq(50, 300, by=30)


## Run the changepoint algorithm on the list of penalty parameters
PELT.result <- vector(mode="list", length=length(pen.list))

# Model: y_t = asym + alpha0*exp(-exp(lgamma)*t) + e_t
for (i in 1:length(pen.list)) {
  # choose a penalty
  pen.i <- pen.list[i]
  
  sm.pelt <- spike.PELT.msl(data=sm.dat, xreg=0, ini.par=NULL, ini.asym=0.01,
                            pen=pen.i, minsl=24, costfun=spike.exp.nlm1,
                            upper.par=c(0.4, 0.5, 1), thresh=0.001, nprune=FALSE)
  
  PELT.result[[i]] <- sm.pelt
  print(i)
}



## ------------------------------------------------------
## (I) Implementing the method in Hocking et al (2013)
## ------------------------------------------------------

## (a) Prepare the data (the annotations)
ID0 <- which(rain1 > 1)
dif0 <- which(diff(ID0) >= 24)
rainID.high <- c(ID0[dif0], rev(ID0)[1])

ID0 <- which(rain1 > 0.5)   # 0.8 for srer and serc
dif0 <- which(diff(ID0) >= 24)
rainID.low <- c(ID0[dif0], rev(ID0)[1])

rain.large1 <- rain.large2 <- rep(0, length(rain1))
rain.large1[rainID.high] <- 1
rain.large2[rainID.low] <- 1

# divide the time axis into regions
wdl <- 24 * 10   # window length
sum(is.na(rain1))

# (1) if there is no missing
region.s <- seq(1, length(rain1), by=wdl)
region.e <- c(region.s[-1]-1, length(rain1))

# (2) if there are missing gaps
# find all the non-missing regions, then break down further if they are too long...
na.rle <- rle(is.na(rain1))
cbind(na.rle$values, na.rle$lengths, cumsum(na.rle$lengths))

re1 <- cumsum(na.rle$lengths)[which(na.rle$values == 0)]
rs1 <- re1 - na.rle$lengths[which(na.rle$values == 0)] + 1

# then break by the windows
wdl0 <- re1 - rs1

wID <- which(wdl0 > wdl)
rs2 <- re2 <- numeric(0)
for (i in wID) {
  rs0 <- seq(rs1[i], re1[i], by=wdl)
  re0 <- c(rs0[-1]-1, re1[i])
  rs2 <- c(rs2, rs0)
  re2 <- c(re2, re0)
}

region.s <- sort(unique(c(rs1, rs2)))  # start ID of the region, 
region.e <- sort(unique(c(re1, re2)))  # end ID of the region

# the annotation
annot <- list()
for (i in 1:length(region.s)) {
  annot1 <- sum(rain.large1[region.s[i]:region.e[i]])
  annot2 <- sum(rain.large2[region.s[i]:region.e[i]])
  rg <- sort(c(annot1, annot2))
  annot[[i]] <- seq(rg[1], rg[2], by=1)
}


## (b) Calculate the loss function
loss.fun <- numeric(0)
for (j in 1:length(PELT.result)) {
  cpts.vec <- rep(0, length(rain1))
  cpts.vec[PELT.result[[j]]$cpt[-length(PELT.result[[j]]$cpt)]] <- 1
  loss.j <- 0
  for (i in 1:length(annot)) {
    ncpt <- sum(cpts.vec[region.s[i]:region.e[i]])
    if (ncpt %in% annot[[i]]) {
      loss.j <- loss.j + 0
    } else {
      loss.j <- loss.j + 1
    }
  }
  # take note
  loss.fun <- c(loss.fun, loss.j)
}



## ------------------------------------------------------
## (II) Implementing the method in Truong et al (2017)
## ------------------------------------------------------

## A function for calculating the cost using labels
spike.exp.label = function(y, ini.par, ini.asym, low.alpha0, upper.par, thresh){
  nt <- length(y)
  segsm <- data.frame(t=1:nt, sm=y)
  
  # the asymptotic parameter
  if (is.null(ini.asym)) ini.asym <- round(0.8 * min(segsm$sm, na.rm=TRUE), 3)
  # the jump and drying parameters
  if (is.null(ini.par)) {
    ini.alpha0 <- round(max(segsm$sm, na.rm=TRUE), 3)
    ini.R <- try(ar(segsm$sm, order.max=1, na.action=na.pass)$ar, silent=TRUE)
    if (length(ini.R) == 1) {
      if ((ini.R < 1) & (ini.R > 0)) {
        ini.lgamma <- round(log(-log(ini.R)), 3)
      } else {
        ini.lgamma <- -5.3  # this is roughly log(-log(0.995))
      }
    } else if (length(ini.R) == 0) {
      ini.lgamma <- -5.3
    }
  } else {
    ini.alpha0 <- ini.par[1]
    ini.lgamma <- ini.par[2]
  }
  
  # fit the exponential model using nlfb
  # the upper limit of the parameter is inherited from the wrapper
  low.alpha <- low.alpha0 + thresh
  fit.nls <- try(
    nlfb(start=list(asym=ini.asym, alpha0=ini.alpha0, lgamma=ini.lgamma), 
         resfn=exp.res1, jacfn=exp.jac1, trace=FALSE, 
         lower=c(0, low.alpha, -20), upper=upper.par,
         data=segsm),
    silent = TRUE
  )
  if (is.list(fit.nls)) {
    coef.vec <- coefficients(fit.nls)
    names(coef.vec) <- c("asym", "alpha0", "lgamma")
    
    # take note of the likelihood regardless
    neglike <- nrow(segsm) * (log(mean(fit.nls$resid^2)) + 1)
    rss <- fit.nls$ssquares
    
  } else {
    # fit a linear model using lgamma
    segsm$expt <- exp(-exp(ini.lgamma) * segsm$t)
    fit.lm <- lm(sm ~ expt, data=segsm)
    
    coef.vec <- coefficients(fit.lm)
    neglike <- nrow(segsm) * (log(mean(fit.lm$residuals^2)) + 1)
    rss <- sum(fit.lm$residuals^2)
  }
  
  output <- list(neglike=neglike, coef.vec=coef.vec, rss=rss)
  return(output)
}


## Using the estimated parameters
spike.exp.cpt <- function(sID, eID, cpts, data, lastchangecoef, type=1){
  
  # find cps for this segment (may be before sID if there's missing)
  ID1 <- which(cpts < sID)
  if (length(ID1 > 0)) {
    cps <- cpts[rev(ID1)[1]] + 1
  } else { cps = 1 }
  
  ID2 <- which(cpts >= eID)
  if (length(ID2 > 0)) {
    cpt <- cpts[ID2[1]]
  } else { cpt = length(data) }
  
  y <- data[cps:cpt]
  nt <- length(y)
  segsm <- data.frame(t=1:nt, sm=y)
  
  # the parameters
  coef.vec <- lastchangecoef[[cpt+1]]
  asym <- coef.vec[1]
  alpha0 <- coef.vec[2]
  lgamma <- coef.vec[3]
  
  # the fitted time series
  t <- segsm$t
  if (type == 1) {
    yhat <- asym + (alpha0 - asym) * exp(-exp(lgamma)*t)
  } else if (type == 2) {
    yhat <- asym + alpha0 * exp(-exp(lgamma)*t)
  }
  sm.resi <- segsm$sm - yhat
  
  ys <- max(c(cps, sID))
  ye <- min(c(cpt, eID))
  subID <- (cps:cpt) %in% (ys:ye)
  neglike <- (ye-ys+1) * (log(mean(sm.resi[subID]^2)) + 1)
  sm.yhat <- yhat[subID]
  
  return(list(neglike=neglike, yhat=sm.yhat))
}


## (a) The overall cost using expert label
# prepare the data (rainfall labels)
ID0 <- which(rain1 >= 1)
dif0 <- which(diff(ID0) >= 24)
rainID <- c(ID0[dif0], rev(ID0)[1])

rain.large <- rep(0, length(rain1))
rain.large[rainID] <- rain1[rainID]

label1 <- which(rain.large > 1)
start.lab0 <- c(1, label1+1)
end.lab0 <- c(label1, length(sm.dat))

# remove the missing gaps
na.rle <- rle(is.na(rain1))
cbind(na.rle$values, na.rle$lengths, cumsum(na.rle$lengths))

# non-missing regions
re1 <- cumsum(na.rle$lengths)[which(na.rle$values == 0)]
rs1 <- re1 - na.rle$lengths[which(na.rle$values == 0)] + 1

# missing regions
re2 <- cumsum(na.rle$lengths)[which(na.rle$values == 1)]
rs2 <- re2 - na.rle$lengths[which(na.rle$values == 1)] + 1

# finding the missing gaps
rm1.list <- rm2.list <- numeric(0)
for (i in 1:length(rs2)) {
  sID <- which(start.lab0 > rs2[i])
  if (length(sID) > 0) {
    rm1 = sID[1]-1
  } else { rm1 = length(start.lab0) }
  
  eID <- which(end.lab0 < rs2[i])
  if (length(eID) > 0) {
    rm2 = rev(eID)[1]+1
  } else { rm2 = 1}
  
  rm1.list <- c(rm1.list, rm1)
  rm2.list <- c(rm2.list, rm2)
}

swap.lab <- numeric(0)
rm.list <- numeric(0)
for (i in 1:length(rm1.list)) {
  t1 <- rm1.list[i]
  t2 <- rm2.list[i]
  
  if (t1 == t2) {
    copy.start <- start.lab0[t1] 
    copy.end <- end.lab0[t1]
    split1 <- c(copy.start, rs2[i])
    split2 <- c(re2[i], copy.end)
    # this is the one to replace
    swap.lab <- rbind(swap.lab, split1, split2)
    rm.list <- c(rm.list, t1)
  } else {
    copy.start <- start.lab0[t1] 
    copy.end <- end.lab0[t2]
    split1 <- c(copy.start, rs2[i])
    split2 <- c(re2[i], copy.end)
    # this is the one to replace
    swap.lab <- rbind(swap.lab, split1, split2)
    rm.list <- c(rm.list, t1:t2)
  }
}

# remove the ones with missing
start.lab <- start.lab0[-rm.list]
end.lab <- end.lab0[-rm.list]
start.lab <- sort(as.vector(c(start.lab, swap.lab[, 1])))
end.lab <- sort(as.vector(c(end.lab, swap.lab[, 2])))

# calculate the overall cost
# initial values
low.alpha0 = 0
ini.par = c(round(max(sm.dat, na.rm=TRUE), 3), -5.3)   
ini.asym = round(0.8 * min(sm.dat, na.rm=TRUE), 3)     

spike.cost <- numeric(0)
for (i in 1:length(start.lab)) {
  s0 <- start.lab[i]
  s1 <- end.lab[i]
  spike.fit <- spike.exp.label(y=sm.dat[s0:s1], ini.par, ini.asym=0.01, 
                               low.alpha0, upper.par=c(0.4, 0.5, 1), 
                               thresh=0.001)
  cost.temp <- spike.fit$neglike  
  # what if it doesn't converge?
  spike.cost <- c(spike.cost, cost.temp)
}

cost.label <- sum(spike.cost)
pcost.label <- cost.label + pen.list * length(start.lab)


## (b) The cost with detected changepoints and estimated parameters
# retrieve the overall cost, excluding the missing gaps
pcost.pelt <- numeric(0)
for (j in 1:length(pen.list)) {
  cpts <- PELT.result[[j]]$cpt
  
  label.cpt <- sm.pelt$cpt[-length(sm.pelt$cpt)]
  
  start.lab0 <- c(1, label.cpt+1)
  end.lab0 <- c(label.cpt, length(sm.dat))
  
  rm1.list <- rm2.list <- numeric(0)
  for (i in 1:length(rs2)) {
    sID <- which(start.lab0 > rs2[i])
    if (length(sID) > 0) {
      rm1 = sID[1]-1
    } else { rm1 = length(start.lab0) }
    
    eID <- which(end.lab0 < rs2[i])
    if (length(eID) > 0) {
      rm2 = rev(eID)[1]+1
    } else { rm2 = 1}
    
    rm1.list <- c(rm1.list, rm1)
    rm2.list <- c(rm2.list, rm2)
  }
  
  swap.lab <- numeric(0)
  rm.list <- numeric(0)
  for (i in 1:length(rm1.list)) {
    t1 <- rm1.list[i]
    t2 <- rm2.list[i]
    
    if (t1 == t2) {
      copy.start <- start.lab0[t1] 
      copy.end <- end.lab0[t1]
      split1 <- c(copy.start, rs2[i])
      split2 <- c(re2[i], copy.end)
      # this is the one to replace
      swap.lab <- rbind(swap.lab, split1, split2)
      rm.list <- c(rm.list, t1)
    } else {
      copy.start <- start.lab0[t1] 
      copy.end <- end.lab0[t2]
      split1 <- c(copy.start, rs2[i])
      split2 <- c(re2[i], copy.end)
      # this is the one to replace
      swap.lab <- rbind(swap.lab, split1, split2)
      rm.list <- c(rm.list, t1:t2)
    }
  }
  
  start.lab <- start.lab0[-rm.list]
  end.lab <- end.lab0[-rm.list]
  start.lab <- sort(as.vector(c(start.lab, swap.lab[, 1])))
  end.lab <- sort(as.vector(c(end.lab, swap.lab[, 2])))
  
  spike.cost.cpt <- numeric(0)
  for (i in 1:length(start.lab)) {
    s0 <- start.lab[i]
    s1 <- end.lab[i]
    spike.fit <- spike.exp.cpt(y=sm.dat[s0:s1], ini.par, ini.asym=0.01, 
                               low.alpha0, upper.par=c(0.4, 0.5, 1), 
                               thresh=0.001)
    cost.temp <- spike.fit$neglike  
    spike.cost.cpt <- c(spike.cost.cpt, cost.temp)
  }
  
  pcost.j <- sum(spike.cost.cpt) + pen.list[j] * length(start.lab)
  pcost.pelt <- c(pcost.pelt, pcost.j)
}


## (c) excess penalised risk
excess.risk = pcost.pelt - pcost.label

