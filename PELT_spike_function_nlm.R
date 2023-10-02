## -----------------------------------------
## PELT functions using the nlmrt package
## -----------------------------------------

## Estimating the changepoints and exponential decay model parameters
## using PELT algorithm
## Solving non-linear least square using LM algorithm (Nash variant)
## User defined upper limits

## Two versions of the exponential decay model
## version 1: y_t = asym + (alpha0 - asym)*exp(-exp(lgamma)*t) + e_t
## version 2: y_t = asym + alpha0*exp(-exp(lgamma)*t) + e_t


library(nlmrt)


## The residual and Jacobian functions
exp.res1 <- function(par, data) {
  # data = data.frame(sm, t, ...)
  # par = c(asym, alpha0, lgamma)
  res <- par[1] + (par[2] - par[1]) * exp(-exp(par[3])*data$t) - data$sm
  return(res)
}

exp.jac1 <- function(par, data) {
  # the Jacobian evaluated at each data point
  Jac <- matrix(0, nrow=nrow(data), ncol=length(par))
  eexp.gamma <- exp(-exp(par[3])*data$t)
  exp.gamma <- -exp(par[3])*data$t
  Jac[, 1] <- 1 - eexp.gamma
  Jac[, 2] <- eexp.gamma
  Jac[, 3] <- (par[2] - par[1]) * eexp.gamma * (exp.gamma)
  return(Jac)
}

exp.res2 <- function(par, data) {
  # data = data.frame(sm, t, ...)
  # par = c(asym, alpha0, lgamma)
  res <- par[1] + par[2] * exp(-exp(par[3])*data$t) - data$sm
  return(res)
}

exp.jac2 <- function(par, data) {
  # the Jacobian evaluated at each data point
  Jac <- matrix(0, nrow=nrow(data), ncol=length(par))
  eexp.gamma <- exp(-exp(par[3])*data$t)
  exp.gamma <- -exp(par[3])*data$t
  Jac[, 1] <- 1
  Jac[, 2] <- eexp.gamma
  Jac[, 3] <- par[2] * eexp.gamma * exp.gamma
  return(Jac)
}


## The cost function
spike.exp.nlm1 = function(y, ini.par, ini.asym, low.alpha0, upper.par, thresh){
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
    if (coef.vec[2] > coef.vec[1]) { 
      sm.hat <- coef.vec[1] + (coef.vec[2] - coef.vec[1]) * exp(-exp(coef.vec[3])*segsm$t)
      neglike <- nrow(segsm) * (log(mean(fit.nls$resid^2)) + 1)
      rss <- fit.nls$ssquares
    } else {
      sm.hat <- NA
      neglike <- 1e20
      rss <- fit.nls$ssquares
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


spike.exp.nlm2 = function(y, ini.par, ini.asym, low.alpha0, upper.par, thresh){ 
  nt <- length(y)
  segsm <- data.frame(t=1:nt, sm=y)
  
  # the asymptotic parameter
  if (is.null(ini.asym)) ini.asym <- round(0.8 * min(segsm$sm, na.rm=TRUE), 3)
  # the jump and drying parameters
  if (is.null(ini.par)) {
    ini.alpha0 <- round(max(segsm$sm, na.rm=TRUE), 3) - ini.asym
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
  fit.nls <- try(
    nlfb(start=list(asym=ini.asym, alpha0=ini.alpha0, lgamma=ini.lgamma), 
         resfn=exp.res2, jacfn=exp.jac2, trace=FALSE, 
         lower=c(0, thresh, -20), upper=upper.par,
         data=segsm),
    silent = TRUE
  )
  if (is.list(fit.nls)) {
    coef.vec <- coefficients(fit.nls)
    names(coef.vec) <- c("asym", "alpha0", "lgamma")
    if (coef.vec[1] + coef.vec[2] > low.alpha0 + thresh) {
      sm.hat <- coef.vec[1] + coef.vec[2] * exp(-exp(coef.vec[3])*segsm$t)
      neglike <- nrow(segsm) * (log(mean(fit.nls$resid^2)) + 1)
      rss <- fit.nls$ssquares
    } else {
      sm.hat <- NA
      neglike <- 1e20
      rss <- fit.nls$ssquares
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


## The PELT iteration
## To initialise, input data, initial parameter values (ini.par, ini.asym), 
## penalty parameter of PELT (pen), minimum segment length (minsl),
## the type of cost function (spike.exp.nlm1 or spike.exp.nlm2), 
## upper limits of parameters and threshold of the jumpsize
spike.PELT.msl <- function(data, xreg=0, ini.par, ini.asym, pen, minsl,
                           costfun, upper.par=NULL, thresh=0, nprune=FALSE){
  # get the cost function
  spike.exp <- costfun
  
  # set the upper limit
  if (is.null(upper.par)) upper.par <- c(0.5, 0.75, 1)
  
  # storage
  n <- length(data)
  lastchangecpts <- rep(NA, n)
  lastchangelike <- c(-pen, rep(NA, n))
  lastchangecoef <- vector(mode="list", length=n+1)
  lastchangecoef[[1]] <- NULL
  lastchangefit <- c(0, rep(NA, n))
  
  # the first few time points
  low.alpha0 = 0  
  for (tstar in minsl:(2*minsl-1)) {
    spike.fit <- spike.exp(y=data[1:tstar], ini.par, ini.asym, low.alpha0, upper.par, thresh)
    lastchangelike[tstar+1] <- spike.fit$neglike
    lastchangecpts[tstar] <- 0
    coef.vec <- c(spike.fit$coef.vec, low.alpha0)
    if (is.null(names(spike.fit$coef.vec))) {names(coef.vec) <- c("null", "low.alpha0")} else
    {names(coef.vec) <- c(names(spike.fit$coef.vec), "low.alpha0")}
    lastchangecoef[[tstar+1]] <- coef.vec
    lastchangefit[tstar+1] <- spike.fit$last.hat
  }
  
  # the rest of the time iterations
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
      spike.fit <- spike.exp(y=data[(checklist[i]+1):tstar], ini.par, ini.asym, 
                             low.alpha0, upper.par, thresh)
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
            # if K is not 0, use tmplike[j] + K <= lastchangelike[tstar+1]+pen
            checklist.remove[j] <- n+2  # keep
          }
          if (tmplike[j] > lastchangelike[tstar+1]+pen) {
            # delay pruning
            # if (checklist.remove[j] > n+1) checklist.remove[j] <- tstar+minsl-1
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
  
  # return the list of object
  PELT.obj <- list(cpt=cpt, lastchangecoef=lastchangecoef,
                   lastchangefit=lastchangefit, 
                   lastchangelike=lastchangelike)
  return(PELT.obj)
}


## Reconstruct the fitted ts from the coefficient list
spike.PELT.yhat <- function(cpt, data, xreg=0, lastchangecoef, type=1){
  
  sm.hat <- rep(NA, times=length(data))
  sm.resi <- rep(NA, times=length(data))
  
  for (i in 1:(length(cpt)-1)) {
    sID <- cpt[i] + 1
    eID <- cpt[i+1]
    y <- data[sID:eID]
    nt <- length(y)
    segsm <- data.frame(t=1:nt, sm=y)
    
    # the parameters
    coef.vec <- lastchangecoef[[eID+1]]
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
    sm.hat[sID:eID] <- yhat
    sm.resi[sID:eID] <- segsm$sm - yhat
  }
  
  PELT.yhat <- list(sm.hat=sm.hat, sm.resi=sm.resi)
  return(PELT.yhat)
}


