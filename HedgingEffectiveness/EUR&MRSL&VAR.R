
rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rio)

data <- import("WeatherIndices&Yields.xlsx")
  
############ Weather indices and yield loss 

A <- as.matrix(data[,4:25])

class(A) <- "numeric"

## whole-cycle
gdd0 <- A[, 8]
cri0 <- A[, 9]
#rdi0 <- A[, 10]

## Stage 1
gdd1 <- A[, 11]
cri1 <- A[, 12]
#rdi1 <- A[, 13]

# Stage 2
gdd2 <- A[, 14]
cri2 <- A[, 15]
#rdi2 <- A[, 16]

# Stage 3 
gdd3 <- A[, 17]
cri3 <- A[, 18]
#rdi3 <- A[, 19]

# Stage 4 
gdd4 <- A[, 20]
cri4 <- A[, 21]
#rdi4 <- A[, 22]

yieldloss <- A[, 7]

## data scaling
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

x01 <- range01(gdd0)
x02 <- range01(cri0)
# x02 <- range01(rdi0)

x11 <- range01(gdd1)
x12 <- range01(cri1)
# x12 <- range01(rdi1)

x21 <- range01(gdd2)
x22 <- range01(cri2)
# x22 <- range01(rdi2)

x31 <- range01(gdd3)
x32 <- range01(cri3)
# x32 <- range01(rdi3)

x41 <- range01(gdd4)
x42 <- range01(cri4)
# x42 <- range01(rdi4)


# exponential utility
ALPHA <- c(0.0052, 0.008, 0.0103)
alpha <- ALPHA[1]
ylexp <- exp(alpha*yieldloss)
ye <- range01(ylexp)

source("basic-functions.R")

source("Train-P.R")

source("Train-W.R")

################# whole-cycle #################

## nseg is the number of knots 
## risk aversion parameter 0.0052, nseg = 21
## risk aversion parameter 0.008 and 0.0103
## nseg = 30

PSAW <- WCPSA(ye, x01, x02, 21)

Fitw <- as.numeric(PSAW$Fit)

Fitw_back <- back(Fitw, ylexp)
Fitw_exp <- (1/alpha)*log(Fitw_back)

Exy.fun.w = Fitw_exp

Mw <- PSAW$M
bw <- PSAW$b
Xbetaw <- Mw[,1:4]%*%bw[1:4]
Xbetaw_back <- back(Xbetaw, ylexp)

#################### phase-division ###################

PSAP <- PDPSA(ye, x11, x12, x21, x22, x31, x32, x41, x42, 39)

Fitp <- PSAP$Fit

Fitp_back <- back(Fitp, ylexp)
Fitp_exp <- (1/alpha)*log(Fitp_back)

Exy.fun.p = Fitp_exp

Mp <- PSA$M
bp <- PSA$b
Xbetap <- Mp[,1:16]%*%bp[1:16]
Xbetap_back <- back(Xbetap, ylexp)

## EU Ratio

XbetaR <- Xbetap_back/Xbetaw_back


## Optimal Indemnity

Istarfun = function(Exy.fun, eta, M){
  pmin(pmax(Exy.fun + eta, 0), M)
}

objFun = function(Exy.fun, eta, M, P){
  res = abs(mean(Istarfun(Exy.fun, eta, M)) - P)
  res
}

## given a combination of P and M

big.n = 0
eta.seq = seq(-10, 10, by=0.1)
for(eta.t in eta.seq)
  big.n = c(big.n, objFun(Exy.fun.w, eta.t, M = 35, P = 9))
big.n = big.n[-1]
etastar.w <- optim(par = 1, objFun, Exy.fun = Exy.fun.w, M = 35, P = 9)$par
Istarmat.w <- Istarfun(Exy.fun.w, eta = etastar.w, M = 35)

big.n = 0
eta.seq = seq(-10, 10, by=0.1)
for(eta.t in eta.seq)
big.n = c(big.n, objFun(Exy.fun.p, eta.t, M = 35, P = 9))
big.n = big.n[-1]
etastar.p <- optim(par = 1, objFun, Exy.fun = Exy.fun.p, M = 35, P = 9)$par
Istarmat.p <- Istarfun(Exy.fun.p, eta = etastar.p, M = 35)

## MRSL mean root square loss
YL <- A[, 7]

YLC <- matrix(NA, nrow = 23, ncol = 96)
YLCM <- vector(mode="logical", length=96)

for (i in (1:96)){
  a = (i-1)*23+1
  b = 23*i
  YLC[,i] <- YL[a:b]
}

for (jj in (1:96)){
  YLCM[jj] <- mean(YLC[,jj])
}

## whole-cycle
MRSLO_w <- MRSLW_w <- vector(mode="logical", length=96)
ISTARW <- matrix(NA, nrow = 23, ncol = 96)

for (i in (1:96)){
  a = (i-1)*23+1
  b = 23*i
  ISTARW[,i] <- Istarmat.w[a:b]
}

for (i in 1:96){
MRSLO_w[i] <- mean(max(0,(YLCM[i] - YLC[,i]))^2)
## VAR, variability of Revenue, is defined as mean((YLCM[i] - YLC[,i])^2)
MRSLW_w[i] <- mean(max(0,(YLCM[i] - YLC[,i] - ISTARW[,i] + 9))^2)
}

PC_w <- (sqrt(MRSLW_w) - sqrt(MRSLO_w))/sqrt(MRSLO_w)

## phase-division
MRSLO_p <- MRSLW_p <- vector(mode="logical", length=96)
ISTARP <- matrix(NA, nrow = 23, ncol = 96)

for (i in (1:96)){
  a = (i-1)*23+1
  b = 23*i
  ISTARP[,i] <- Istarmat.p[a:b]
}

for (i in 1:96){
  MRSLO_p[i] <- mean(max(0,(YLCM[i] - YLC[,i]))^2)
  MRSLW_p[i] <- mean(max(0,(YLCM[i] - YLC[,i] - ISTARP[,i] + 9))^2)
}

PC_p <- (sqrt(MRSLW_p) - sqrt(MRSLO_p))/sqrt(MRSLO_p)

PC <- cbind(as.vector(PC_w), as.vector(PC_p))

library(writexl)
PC <- as.data.frame(PC)
write_xlsx(PC, "C:\\...\\filename.xlsx")

## the ratio between two expectations

fenzi <- exp((-1)*alpha*Istarmat.p)

fenmu <- exp((-1)*alpha*Istarmat.w)

EURatio <- XbetaR*(fenzi/fenmu)

for (i in (1:96)){
  a = (i-1)*23+1
  b = 23*i
EuRatiomean[,i] <- EURatio[a:b]
}

Eumean_county <- vector(mode="logical", length=96)
for (jj in (1:96)){
  Eumean_county[jj] <- mean(EuRatiomean[,jj])
}

# check alpha to name the EuRatiomean
alpha

EuRatiomean0103 <- as.vector(Eumean_county)

## after obtain three EuRatiomean, 
## combine them together
EuRatiomean_county <- cbind(EuRatiomean0052, EuRatiomean008, EuRatiomean0103)

library(writexl)
EuRatiomean_county <- as.data.frame(EuRatiomean_county)
write_xlsx(EuRatiomean_county, "C:\\...\\filename.xlsx")

############ across counties
library(rio)
PC <- import("EuRatiomean.xlsx")
PC <- as.matrix(PC)
class(PC) <- "numeric"

head(PC)

Region <- PC[,3]
mean0052 <- PC[,4]
mean008 <- PC[,5]
mean0103 <- PC[,6]

W0052 <- PC[,4]
P0052<- PC[,5]
W008 <- PC[,6]
P008 <- PC[,7]
W0103 <- PC[,8]
P0103 <- PC[,9]

mean(P0052)
mean(P008)
mean(P0103)

colorgroup <- c(rgb(26/255,51/255,153/255), rgb(49/255,98/255,173/255), rgb(75/255,152/255,196/255), rgb(125/255,203/255,212/255), rgb(195/255,230/255,196/255), rgb(219/255,216/255,134/255), rgb(194/255,162/255,66/255), rgb(163/255,98/255,28/255), rgb(128/255,25/255,0))

library(grid)
library(gridExtra)
library(scales)
library(memisc)

pdf(file="C:\\EuRatio_plot.pdf", height = 3, width=8)
par(mfrow = c(1,3))
plot(Region[1:96], P0052[1:96], yaxt="n", pch = 21, xlab="Agricultural District", ylab="MRSL Percentage Change", cex.lab = 1, cex.axis = 1, ylim = c(-0.7, 0.3), cex = 1.5, bg=colorgroup[as.factor(Region)])
axis(2, at=pretty_breaks(n=5)(P0052[1:96]), lab=paste0(pretty(P0052[1:96]) * 100, " %"), las=TRUE)
abline(h=1, col="black", lty=2)
title(main=expression(paste(alpha, " = ", 0.0052, "_P")), cex.main = 1.2, adj = 0)

plot(Region[1:96], P008[1:96], yaxt="n", pch = 21, xlab="Agricultural District", ylab = "MRSL Percentage Change", cex.lab = 1, cex.axis = 1, ylim = c(-0.7, 0.3), cex = 1.5, bg=colorgroup[as.factor(Region)])
axis(2, at=pretty_breaks(n=5)(P008[1:96]), lab=paste0(pretty(P008[1:96]) * 100, " %"), las=TRUE)
abline(h=1, col="black", lty=2)
title(main=expression(paste(alpha, " = ", 0.008, "_P")), cex.main = 1.2, adj = 0)

plot(Region[1:96], P0103[1:96], yaxt="n", pch = 21, xlab="Agricultural District", ylab = "MRSL Percentage Change", cex.lab = 1, cex.axis = 1, cex = 1.5, ylim = c(-0.7, 0.3), bg=colorgroup[as.factor(Region)])
axis(2, at=pretty_breaks(n=5)(P0103[1:96]), lab=paste0(pretty(P0103[1:96]) * 100, " %"), las=TRUE)
abline(h=1, col="black", lty=2)
title(main=expression(paste(alpha, " = ", 0.0103, "_P")), cex.main = 1.2, adj = 0)
dev.off()

