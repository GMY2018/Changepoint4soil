## -------------------------------------------
## NEON soil moisture time series analysis
## -------------------------------------------

## An example of adding covariates

source("PELT_soil_function_nlm.R")
load("Rain_SM_data_SRER.RData")


tID <- 35000:42000
plot(neon.sm.top$time[tID], neon.sm.top$hloc_2[tID], type="l", 
     xlab="time", ylab="soil moisture")

## Covariate: rainfall
rain <- neon.rain$precip
rain[is.na(rain)] <- 0

# only keep the significant precipitation event
quantile(rain, probs=c(0.95, 0.99, 0.995, 0.996, 0.997, 0.998, 0.999))
rain2 <- rain
rain2[rain <= 0.95] <- 0      

diff.t <- diff(which(rain2 > 0))
rain.t <- cumsum(diff.t)[which(diff.t > 10)]
rain.dat0 <- rep(0, length(rain2))
rain.dat0[rain.t] <- 1
rain.dat <- rain.dat0[tID]

smx <- rain.dat[seq(1, length(sm.dat), by=2)]
smy <- sm.dat[seq(1, length(sm.dat), by=2)]

par(mfrow=c(2, 1))
plot(smy, type="l")
plot(smx, type="h")



## New cost functions for using covariates
## Version 1: including rainfall data as indicator functions
spike.exp.xreg1 = function(y, x, ini.par, ini.asym, low.alpha0, 
                           upper.par=NULL, thresh=0){ 
  # the data frame
  nt <- length(y)
  brk <- which(x > 0)
  brk1 <- c(brk[-1]-1, nt)
  if (length(brk) > 0) {
    rain.bin <- matrix(0, nrow=nt, ncol=length(brk))
    for (i in 1:ncol(rain.bin)) {
      rain.bin[brk[i]:nt, i] <- 1  
    }
    segsm <- data.frame(t=1:nt, sm=y, rain.bin)
    names(segsm) <- c("t", "sm", paste0("rain.b", 1:ncol(rain.bin)))
  } else {
    segsm <- data.frame(t=1:nt, sm=y)
  }
  
  if (is.null(ini.par)) {
    ini.alpha0 <-round(max(segsm$sm, na.rm=TRUE), 3)
    ini.R <- try(ar(segsm$sm, order.max=1, na.action=na.pass)$ar, silent=TRUE)
    if (length(ini.R) == 1) {
      if ((ini.R < 1) & (ini.R > 0)) {
        ini.lgamma <- round(log(-log(ini.R)), 3)
      } else {
        ini.lgamma <- -5.3     # this is roughly log(-log(0.995))
      }
    } else if (length(ini.R) == 0) {
      ini.lgamma <- -5.3
    }
  } else {
    ini.alpha0 <- ini.par[1]
    ini.lgamma <- ini.par[2]
  }
  # the asymptotic parameter
  if (is.null(ini.asym)) ini.asym <- round(0.8 * min(segsm$sm, na.rm=TRUE), 3)
  
  # use nls() to fit the model
  fit.nls <- spike.nls(segsm, ini.asym, ini.alpha0, ini.lgamma, low.alpha0)
  if (is.list(fit.nls)) {
    coef.vec <- coefficients(fit.nls)
    if (coef.vec[1] + coef.vec[2] > low.alpha0 + thresh) {
      sm.hat <- coef.vec[1] + coef.vec[2] * exp(-exp(coef.vec[3])*segsm$t)
      # neglike <- nrow(segsm) * (log(mean(residuals(fit.nls)^2, na.rm=TRUE)) + 1)
      neglike <- nrow(segsm) * (log(mean((sm.hat - segsm$sm)^2, na.rm=TRUE)) + 1)
      rss <- mean(residuals(fit.nls)^2)
    } else {
      sm.hat <- NA
      neglike <- 1e20
      rss <- NA
    }
  } else {
    coef.vec <- NA
    sm.hat <- NA
    neglike <- 2e20
    rss <- NA
  }
  
  # return(neglike)
  output <- list(neglike=neglike, coef.vec=coef.vec, rss=rss, last.hat=rev(sm.hat)[1])
  return(output)
}


## Version 2: including rainfall data as step function
spike.exp.xreg2 = function(y, x, ini.par, ini.asym, low.alpha0, 
                           upper.par=NULL, thresh=0){ 
  # the data frame
  nt <- length(y)
  brk <- which(x > 0)
  brk1 <- c(brk[-1]-1, nt)
  if (length(brk) > 0) {
    rain.binx <-  rep(0:(length(brk)-1), times=diff(c(0, brk1)))
    segsm <- data.frame(t=1:nt, sm=y, rain.binx)
    names(segsm) <- c("t", "sm", "rain")
  } else {
    segsm <- data.frame(t=1:nt, sm=y)
  }
  
  if (is.null(ini.par)) {
    ini.alpha0 <-round(max(segsm$sm, na.rm=TRUE), 3)
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
  # the asymptotic parameter
  if (is.null(ini.asym)) ini.asym <- round(0.8 * min(segsm$sm, na.rm=TRUE), 3)
  
  # use nls() to fit the model
  fit.nls <- spike.nls(segsm, ini.asym, ini.alpha0, ini.lgamma, low.alpha0)
  if (is.list(fit.nls)) {
    coef.vec <- coefficients(fit.nls)
    if (coef.vec[1] + coef.vec[2] > low.alpha0 + thresh) {
      sm.hat <- coef.vec[1] + coef.vec[2] * exp(-exp(coef.vec[3])*segsm$t)
      # neglike <- nrow(segsm) * (log(mean(residuals(fit.nls)^2, na.rm=TRUE)) + 1)
      neglike <- nrow(segsm) * (log(mean((sm.hat - segsm$sm)^2, na.rm=TRUE)) + 1)
      rss <- mean(residuals(fit.nls)^2)
    } else {
      sm.hat <- NA
      neglike <- 1e20
      rss <- NA
    }
  } else {
    coef.vec <- NA
    sm.hat <- NA
    neglike <- 2e20
    rss <- NA
  }
  
  output <- list(neglike=neglike, coef.vec=coef.vec, rss=rss, last.hat=rev(sm.hat)[1])
  return(output)
}


## the model fitting function
spike.nls = function(segsm, ini.asym, ini.alpha0, ini.lgamma, low.alpha0) {
  if (ncol(segsm) > 2) {
    segsm$sm2 <- segsm$sm - ini.asym
    x.name <- names(segsm)[-c(1:2, ncol(segsm))]
    x.coef <- paste0("b", 1:length(x.name))
    x.fun <- paste(paste(x.coef, x.name, sep="*"), collapse=" + ")
    f.nls <- as.formula(paste("sm2 ~ asym + (alpha0 +", x.fun, ") * exp(-exp(lgamma)*t)"))
    f.ini <- vector(mode="list", length=length(x.coef))
    f.ini[1:length(x.coef)] <- 0   # initial values of the covariates
    names(f.ini) <- x.coef
    fit.nls <- try(
      nls(f.nls, data=segsm, algorithm="port",
          start=c(asym=ini.asym, alpha0=ini.alpha0, lgamma=ini.lgamma, f.ini),
          lower=c(0, low.alpha0, -1e10, rep(0, length(f.ini)))),
      silent = TRUE
    )
  } else {
    segsm$sm2 <- segsm$sm - ini.asym
    f.nls <- as.formula("sm ~ asym + alpha0 * exp(-exp(lgamma)*t)")
    fit.nls <- try(
      nls(f.nls, data=segsm, algorithm="port",
          start=list(asym=ini.asym, alpha0=ini.alpha0, lgamma=ini.lgamma),
          lower=c(0, low.alpha0, -1e10)),
      silent = TRUE
    )
  }
  return(fit.nls)
}



## PELT iterations (I didn't write a wrapper function for this)
data = smy
xreg = smx      
pen = 250
ini.par = NULL
ini.asym = 0.01
low.alpha0 = 0.01
nprune = FALSE
minsl = 24
thresh = 0.01

spike.exp = spike.exp.xreg2

# storage
n <- length(data)
lastchangecpts <- rep(NA, n)
lastchangelike <- c(-pen, rep(NA, n))
lastchangecoef <- vector(mode="list", length=n+1)
lastchangecoef[[1]] <- NULL
lastchangefit <- c(0, rep(NA, n))

# the first few time points
for (tstar in minsl:(2*minsl-1)) {
  spike.fit <- spike.exp(y=data[1:tstar], x=xreg[1:tstar], ini.par, ini.asym, 
                         low.alpha0, upper.par=NULL, thresh)
  lastchangelike[tstar+1] <- spike.fit$neglike
  lastchangecpts[tstar] <- 0
  coef.vec <- c(spike.fit$coef.vec, low.alpha0)
  if (is.null(names(spike.fit$coef.vec))) {names(coef.vec) <- c("null", "low.alpha0")} else
  {names(coef.vec) <- c(names(spike.fit$coef.vec), "low.alpha0")}
  lastchangecoef[[tstar+1]] <- coef.vec
  lastchangefit[tstar+1] <- spike.fit$last.hat
}
# The first few rows computed first because of the minimum segment length
# The cost is C(y_{1:tstar})

noprune = NULL
checklist = 0
checklist.remove <- n+2

for (tstar in (2*minsl):n) {
  checklist <- c(checklist, tstar-minsl)
  checklist.remove <- c(checklist.remove, n+2)
  nchecklist <- length(checklist)
  
  tmplike <- 1:nchecklist
  tmpcost <- 1:nchecklist    # vector of C(y_{t+1:tstar})
  tmplast <- lastchangelike[checklist+1]   # vector of F(t)
  tmpcoef <- vector(mode="list", length=nchecklist)   
  tmpfit <- 1:nchecklist
  
  for (i in 1:nchecklist) {
    # set the lower limit of alpha (for positive spikes)
    lastfit <- lastchangefit[checklist[i]+1]
    if (!is.na(lastfit)) {
      low.alpha0 <- lastfit
    } else {
      low.alpha0 <- min(data, na.rm=TRUE)
    }
    spike.fit <- spike.exp(y=data[(checklist[i]+1):tstar], 
                           x=xreg[(checklist[i]+1):tstar],
                           ini.par, ini.asym, low.alpha0, upper.par=NULL, thresh)
    tmpcost[i] <- spike.fit$neglike
    tmplike[i] <- tmplast[i] + tmpcost[i] + pen 
    # this is F(t) + C(y_{t+1:tstar}) + beta
    coef.vec <- c(spike.fit$coef.vec, low.alpha0)
    if (is.null(names(spike.fit$coef.vec))) {names(coef.vec) <- c("null", "low.alpha0")} else
    {names(coef.vec) <- c(names(spike.fit$coef.vec), "low.alpha0")}
    tmpcoef[[i]] <- coef.vec
    tmpfit[i] <- spike.fit$last.hat
  }
  
  # upeate the likelihood and the last cpt
  if (any(tmplike < 1e19)) {
    # if there is a t in checklist where F(t) and C(y_{t+1:tstar}) are both finite
    # use the usual procedure
    lastchangelike[tstar+1] <- min(tmplike, na.rm=TRUE)
    lastchangecpts[tstar] <- checklist[which.min(tmplike)]
    for (j in 1:nchecklist) {
      if (tmplike[j] < 1e19) {  
        if (tmplike[j] <= lastchangelike[tstar+1]+pen) {
          checklist.remove[j] <- n+2  # keep
        }
        if (tmplike[j] > lastchangelike[tstar+1]+pen) {
          # delay pruning
          if (checklist.remove[j] > n+1) checklist.remove[j] <- checklist[j] + 2*minsl
        }
      } else {
        if (tmplast[j] > 1e19) {
          checklist.remove[j] <- 0  # prune
        } else {  
          # delay pruning (keep for a bit longer)
          if (checklist.remove[j] > n+1) checklist.remove[j] <- checklist[j] + 2*minsl
        }
      }
    }
    rID <- checklist.remove > tstar
    checklist <- checklist[rID]
    checklist.remove <- checklist.remove[rID]
    # record the coefficients and the last fitted values
    lastchangecoef[[tstar+1]] <- tmpcoef[tmplike == lastchangelike[tstar+1]][[1]]
    lastchangefit[tstar+1] <- tmpfit[tmplike == lastchangelike[tstar+1]][1]
    
  } else if (any(tmplast < 1e19)) {
    # if all tmplike are infinite, but some of the tmplast are finite, then we keep them, 
    # and prune the t that has infinite tmplast, for they should never be the last optimal
    like.vec <- tmplike[tmplast < 1e19]
    cpt.vec <- checklist[tmplast < 1e19]
    lastchangelike[tstar+1] <- min(like.vec, na.rm=TRUE)
    lastchangecpts[tstar] <- cpt.vec[which.min(like.vec)]
    for (j in 1:nchecklist) {
      if (tmplast[j] > 1e19) {
        checklist.remove[j] <- 0
      } else if ((tmplast[j] < 1e19)) {
        # delay pruning  (keep for a bit longer)
        if (checklist.remove[j] > n+1) checklist.remove[j] <- checklist[j] + 2*minsl
      }
    }
    rID <- checklist.remove > tstar
    if (sum(rID) > 0) {
      checklist <- checklist[rID]
      checklist.remove <- checklist.remove[rID]
    } else if (sum(rID) == 0) {
      checklist <- cpt.vec
      checklist.remove <- rep(tstar+1, length(cpt.vec))
    }
    lastchangecoef[[tstar+1]] <- tmpcoef[tmplike == lastchangelike[tstar+1]][[1]]
    lastchangefit[tstar+1] <- tmpfit[tmplike == lastchangelike[tstar+1]][1]
    
  }
  
  if (nprune == TRUE) {
    noprune = c(noprune, length(checklist))
  }
} 

fcpt = NULL
last = n
while (last > 0) {
  fcpt = c(fcpt, last)
  last = lastchangecpts[last]
}
cpt <- rev(fcpt)

