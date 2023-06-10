##############################################
# File: pd-MGCV.r 
# Author: [Jing Zou] 
# Created: [2023.06.10] 
# Description: [estimate phase-division models using mgcv package gam()] 
##############################################

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rio)

data <- import("WeatherIndices&Yields.xlsx")

############ Weather indices and yield loss 

A <- as.matrix(data[,4:25])

class(A) <- "numeric"

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


############ data scaling

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

## under quadratic utility, response is yield loss.
y <- range01(yieldloss)

## under exponential utility, reponse is ylexp
#ALPHA <- c(0.0052, 0.008, 0.0103)
#alpha <- ALPHA[1]
#ylexp <- exp(alpha*yieldloss)
# y <- range01(ylexp)

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

### Benchmark method mgcv::gam ###

library(mgcv)

## quadradti utility case

illisoy <- data.frame(x11, x12, x21, x22, x31, x32, x41, x42, y)

fit.pdm <- gam(y~te(x11, x12, bs = "ps") + te(x21, x22, bs = "ps") + te(x31, x32, bs = "ps") + te(x41, x42, bs = "ps"), data = illisoy, method = "REML")


## adj.R^2

fitted.pdm <- fit.pdm$fitted.values
expvar01 <- sum((y-fitted.pdm)^2)
expvar02 <- sum((y-mean(y))^2)
R2pdm = 1 - expvar01/expvar02
n = length(y)
k=8

adR2pdm = 1- (1-R2pdm)*(n-1)/(n-k-1)
adR2pdm

## under quadratic utility case, response y is the yield loss

# function to convert fitted values back to original scale
back <- function(x, y){
  z = x * (max(y)-min(y)) + min(y)
  return(z)
}

y_back <- back(y, yieldloss)
Fit_pdmv_back <- back(fitted.pdm, yieldloss)
RMSE_pdmv <- sqrt(mean((y_back - Fit_pdmv_back)^2))
RMSE_pdmv

## under exponential utility case, 
## response is exp(alpha*yieldloss)

ALPHA <- c(0.0052, 0.008, 0.0103)
## take the low level risk aversion for example
alpha <- ALPHA[1]
ylexp <- exp(alpha*yieldloss)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
ye <- range01(ylexp)

illisoye <- data.frame(x11, x12, x21, x22, x31, x32, x41, x42, ye)

fit.pdme <- gam(ye~te(x11, x12, bs = "ps") + te(x21, x22, bs = "ps") + te(x31, x32, bs = "ps") + te(x41, x42, bs = "ps"), data = illisoye, method = "REML")
fitted.pdme <- fit.pdme$fitted.values

# RMSE
ye_back <- back(ye,ylexp)
y_exp <- (1/alpha)*log(ye_back)
Fit_pdme_back <- back(fitted.pdme, ylexp)
Fit_pdme <- (1/alpha)*log(Fit_pdme_back)
RMSE_pdme <- sqrt(mean((y_exp - Fit_pdme)^2))
RMSE_pdme

# adj.R^2
fitted.pdme <- fit.pdme$fitted.values
expvar01e <- sum((ye-fitted.pdme)^2)
expvar02 <- sum((ye-mean(ye))^2)
R2pdme = 1 - expvar01e/expvar02
n = length(ye)
k = 8

adR2pdme = 1- (1-R2pdme)*(n-1)/(n-k-1)
adR2pdme
