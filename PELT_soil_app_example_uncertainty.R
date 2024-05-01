## --------------------------------------------------
## Extracting the parameters and their uncertainty
## --------------------------------------------------

library(nlmrt)
library(matrixcalc)

source("PELT_spike_function_nlm.R")

load("DataExample/PELT_SRER_data_loc4.RData")
load("DataExample/PELT_SRER_result_example.RData")


## Read the data
sm.loc4 <- sm.dat1[seq(2, length(time0), by=2)]
st4 <- time0[seq(2, length(time0), by=2)]
rain4 <- rain1[seq(2, length(time0), by=2)] + rain1[seq(1, length(time0)-1, by=2)]


## An example using the changepoint detection analysis 
## on data from NEON site SRER
data = sm.loc4
sm.pelt = PELT.result.loc4[[3]]
st = st4
cpt = c(0, sm.pelt$cpt)
lastchangecoef = sm.pelt$lastchangecoef
lastchangelike = sm.pelt$lastchangelike
thresh = 0.001
ini.asym = 0.01
ini.par = NULL


## Extract the parameters
srer.coef <- sapply(lastchangecoef[cpt+1], FUN=as.vector)
srer.coef <- t(srer.coef)
srer.coef <- cbind(srer.coef[, 1:3], diff(c(0, cpt)))
colnames(srer.coef) <- c("asym", "spike", "lgamma", "length")


## Extracting the uncertainty
## Here we extracting the uncertainty by re-fitting the models and 
## record the standard deviations of the estimates

fmodel <- sm ~ asym + alpha0*exp(-exp(lgamma)*t)
fpar <- c(0.01, 0.1, -5)
names(fpar) <- c("asym", "alpha0", "lgamma")
exp.ss2 <- model2ssfun(fmodel, fpar, funname="exp.ss2")
exp.gr2 <- model2grfun(fmodel, fpar, funname="exp.gr2")
exp.hess2 <- function(par, data){
  Hessian <- array(0, dim=c(length(par), length(par), nrow(data)))
  eexp.gamma <- exp(-exp(par[3])*data$t)
  exp.gamma <- -exp(par[3])*data$t
  Hessian[2, 3, ] <- Hessian[3, 2, ] <- eexp.gamma * exp.gamma
  Hessian[3, 3, ] <- par[2] * eexp.gamma * exp.gamma * (exp.gamma + 1)
  return(-Hessian)
}


Jacobian.list <- vector(mode="list", length=length(cpt)-1)  # for the cross product
Deriv.list <- vector(mode="list", length=length(cpt)-1)
Resid.list <- Resid2.list <- vector(mode="list", length=length(cpt)-1)
sighat <- sighat2 <- rep(NA, length(cpt)-1)
SE.mat <- SE2.mat <- matrix(NA, nrow=length(cpt)-1, ncol=3)
coef.mat <- coef2.mat <- matrix(NA, nrow=length(cpt)-1, ncol=3)

for (i in 1:(length(cpt)-1)) {
  sID <- cpt[i] + 1
  eID <- cpt[i+1]
  y <- data[sID:eID]
  nt <- length(y)
  segsm <- data.frame(t=1:nt, sm=y)
  
  # the parameters
  coef.vec <- lastchangecoef[[eID+1]]
  
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
  
  # re-fit the time series and get the Jacobian
  fit.nls <- try(
    nlfb(start=list(asym=ini.asym, alpha0=ini.alpha0, lgamma=ini.lgamma), 
         resfn=exp.res2, jacfn=exp.jac2, trace=FALSE, 
         lower=c(0, thresh, -20), upper=c(0.4, 0.5, 3),
         data=segsm),
    silent = TRUE
  )
  if (is.list(fit.nls)) {
    Resid.list[[i]] <- fit.nls$resid
    coef.mat[i, ] <- fit.nls$coefficients
    Jacobian.list[[i]] <- crossprod(fit.nls$jacobian)
    nls.sum <- summary(fit.nls)
    SE.mat[i, ] <- nls.sum$SEs
    ssnew <- exp.ss2(fit.nls$coefficients, sm=segsm$sm, t=segsm$t)
    sighat[i] <- ssnew / (nrow(segsm) - 3)
  }
  
  # or using nls() function
  fit.nls2 <- try(
    nls(sm ~ asym + alpha0*exp(-exp(lgamma)*t), data=segsm,
        start=list(asym=ini.asym, alpha0=ini.alpha0, lgamma=ini.lgamma),
        trace=FALSE, algorithm="port",
        lower=c(0, thresh, -20), upper=c(0.4, 0.5, 3)),
    silent = TRUE
  )
  if (is.list(fit.nls2)) {
    Resid2.list[[i]] <- residuals(fit.nls2)
    nls.sum2 <- summary(fit.nls2)
    coef2.mat[i, ] <- nls.sum2$coefficients[, 1]
    SE2.mat[i, ] <- nls.sum2$coefficients[, 2]
    ssnew <- exp.ss2(nls.sum2$coefficients[, 1], sm=segsm$sm, t=segsm$t)
    sighat2[i] <- ssnew / (nrow(segsm) - 3)
  }
}


## Sometimes the fit does not have a converged standard deviation
## Therefore, we combine the result from two different fitting functions to 
## maximise the chance that we have a converged result

diff.mat <- coef.mat - srer.coef[, 1:3]
dID1 <- which(apply(abs(diff.mat) > 0.05, MARGIN=1, FUN=any))
diff.mat <- coef2.mat - srer.coef[, 1:3]
dID2 <- which(apply(abs(diff.mat) > 0.05, MARGIN=1, FUN=any))

stdev.mat <- SE.mat
stdev.mat[dID1, ] <- NA

mID <- which(apply(is.na(stdev.mat), MARGIN=1, FUN=any))
coef.mat[mID, ]    # boundary solutions
cbind(coef.mat[mID, 3], SE2.mat[mID, 3])
rID <- mID[!(mID %in% dID2)]
stdev.mat[rID, ] <- SE2.mat[rID, ]


## Visualise the temporal dynamic of parameters
std.asym <- stdev.mat[,1]
std.spike <- stdev.mat[,2]
std.gamma <- stdev.mat[,3]

asym.vec0 <- srer.coef[, 1]
spike.vec0 <- srer.coef[, 2]
gamma.vec0 <- srer.coef[, 3]


## Using rectangles or polygons to display SEs
col2 <- c("grey60", "grey80")

sm.fit <- spike.PELT.yhat(cpt=cpt, data=data, xreg=0, 
                          lastchangecoef=lastchangecoef, type=2)

pdf(file="Dynamic_parameter_example.pdf", width=9, height=8, pointsize=10)
par(mfrow=c(3, 1), mar=c(4.5, 4.5, 2, 2), cex=1.1)

# the changepoints
plot(st, data, type="l", lwd=2, col="gray30", xlab="Time", ylab="Soil water content",
     ylim=c(0.8*min(data), 1.1*max(data)))
lines(st, sm.fit$sm.hat, col="brown3", lwd=2)
points(st[cpt[-1]], rep(0.9*min(data), times=length(cpt)-1), pch=17, 
       col="gray30", cex=0.6)

# alpha_0, the asymptotic level parameter
std.vec <- rep(std.asym, diff(cpt))
asym.vec <- rep(asym.vec0, diff(cpt))

plot(st, asym.vec, type="n", xlab="Time", ylab="alpha_0", ylim=c(-0.05, 0.15))
yb <- asym.vec0-1.96*std.asym
yb[yb < -0.2] <- -0.2
yt <- asym.vec0+1.96*std.asym
yt[yt > 0.4] <- 0.4
# check if need this
mID <- which(is.na(std.asym))
bID <- which((asym.vec0 <= 0) | (asym.vec0 >= 0.4))
sID <- mID[mID %in% bID]
yt[sID] <- 0.4
yb[sID] <- -0.2
rect(xleft=st[cpt[-length(cpt)]+1], xright=st[cpt[-1]],
     ybottom=yb, ytop=yt, col="grey80", border="grey80", lwd=2)
segments(x0=st[cpt[-length(cpt)]+1], x1=st[cpt[-1]], y0=asym.vec[cpt[-1]],
         col="grey30", lwd=3)


# gamma, the transformed decay parameter
std.vec <- rep(std.gamma, diff(cpt))
gamma.vec <- rep(gamma.vec0, diff(cpt))
tau.vec <- 1/exp(gamma.vec)/24

plot(st, gamma.vec, type="n", xlab="Time", ylab="gamma", ylim=c(-15, 10))
yb <- gamma.vec0-1.96*std.gamma
yb[yb < -25] <- -25
yt <- gamma.vec0+1.96*std.gamma
yt[yt > 10] <- 10
# check if need the following 
mID <- which(is.na(std.gamma))
bID <- which((gamma.vec0 <= -20) | (gamma.vec0 >= 3))
sID <- mID[mID %in% bID]
yt[sID] <- 10
yb[sID] <- -25
rect(xleft=st[cpt[-length(cpt)]+1], xright=st[cpt[-1]],
     ybottom=yb, ytop=yt, col="grey80", border="grey80", lwd=2)
segments(x0=st[cpt[-length(cpt)]+1], x1=st[cpt[-1]], y0=gamma.vec[cpt[-1]],
         col="grey30", lwd=3)

dev.off()


