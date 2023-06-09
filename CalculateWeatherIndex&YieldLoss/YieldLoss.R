
rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rio)

data <- import("WeatherIndices&Yields.xlsx")

############ Weather indices and yield data 

A <- as.matrix(data[,4:25])

class(A) <- "numeric"

## Yield Detrending
library(plm)
Detrend  <- as.data.frame(A)
head(Detrend)
E <- pdata.frame(Detrend, index = 96)
## fixed effect 
fe <- plm(yield ~ year_d, data = E, model = "within")
summary(fe)
r <- as.vector(residuals(fe))

within_intercept(fe)

## take the "overall_intercept" from
## the results of within_intercept(fe)
detyield <- r + overall_intercept

### calculate detrended yield loss separately
library(pracma)
c <- matrix(nrow=23, ncol=96)
d <- matrix(nrow=23, ncol=96)
yieldloss <- vector(mode="logical", length=nrow(A))
for (i in (1:96)){
  a = (i-1)*23+1
  b = 23*i
  c[,i] <- r[a:b]
  d[,i] <- max(c[,i]) - c[,i]
  yieldloss[a:b] <- d[,i]
}

