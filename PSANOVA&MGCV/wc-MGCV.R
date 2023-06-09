
rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rio)

data <- import("WeatherIndices&Yields.xlsx")

############ Weather indices and yield loss 

A <- as.matrix(data[,4:25])

class(A) <- "numeric"

gdd0 <- A[, 8]
cri0 <- A[, 9]
#rdi0 <- A[, 10]

yieldloss <- A[, 7]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

## under quadratic utility, response is yield loss.
y <- range01(yieldloss)

x01 <- range01(gdd0)
x02 <- range01(cri0)
# x02 <- range01(rdi0)

### Benchmark method mgcv::gam ###

library(mgcv)

## quadradti utility case

illisoyw <- data.frame(x01, x02, y)

fit.pdm.w <- gam(y~te(x01, x02, bs = "ps"), data = illisoyw, method = "REML")

## adj.R^2

fitted.pdm.w <- fit.pdm.w$fitted.values
expvar01 <- sum((y-fitted.pdm.w)^2)
expvar02 <- sum((y-mean(y))^2)
R2pdmw = 1 - expvar01/expvar02
n = length(y)
k = 2

adR2pdmw = 1- (1-R2pdmw)*(n-1)/(n-k-1)
adR2pdmw

## under quadratic utility case, response y is the yield loss

# function to convert fitted values back to original scale
back <- function(x, y){
  z = x * (max(y)-min(y)) + min(y)
  return(z)
}

y_back <- back(y, yieldloss)
Fit_wmv_back <- back(fitted.pdm.w, yieldloss)
RMSE_wmv <- sqrt(mean((y_back - Fit_wmv_back)^2))
RMSE_wmv

## under exponential utility case, 
## response is exp(alpha*yieldloss)

ALPHA <- c(0.0052, 0.008, 0.0103)
## take the low level risk aversion for example
alpha <- ALPHA[3]
ylexp <- exp(alpha*yieldloss)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
ye <- range01(ylexp)

illisoywe <- data.frame(x01, x02, ye)

fit.pdmwe <- gam(ye~te(x01, x02, bs = "ps"), data = illisoywe, method = "REML")
fitted.pdmwe <- fit.pdmwe$fitted.values

# RMSE
ye_back <- back(ye,ylexp)
y_exp <- (1/alpha)*log(ye_back)
Fit_pdmwe_back <- back(fitted.pdmwe, ylexp)
Fit_pdmwe <- (1/alpha)*log(Fit_pdmwe_back)
RMSE_pdmwe <- sqrt(mean((y_exp - Fit_pdmwe)^2))
RMSE_pdmwe

# adj.R^2
fitted.pdmwe <- fit.pdmwe$fitted.values
expvar01we <- sum((ye-fitted.pdmwe)^2)
expvar02 <- sum((ye-mean(ye))^2)
R2pdmwe = 1 - expvar01we/expvar02
n = length(ye)
k = 2

adR2pdmwe = 1- (1-R2pdmwe)*(n-1)/(n-k-1)
adR2pdmwe