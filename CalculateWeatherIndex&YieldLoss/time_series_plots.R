rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rio)

data <- import("WeatherIndices&Yields.xlsx")

############ Weather indices and yield data 

A <- as.matrix(data[,4:25])

class(A) <- "numeric"

year <- A[, 2]
## yield data

yield <- A[, 5]
detyield <- A[, 6]
yieldloss <- A[, 7]

### whole-cycle

gdd0 <- A[, 8]
cri0 <- A[, 9]

### phase-division
## Stage 1
gdd1 <- A[, 11]
cri1 <- A[, 12]

# Stage 2
gdd2 <- A[, 14]
cri2 <- A[, 15]

# Stage 3 
gdd3 <- A[, 17]
cri3 <- A[, 18]

# Stage 4 
gdd4 <- A[, 20]
cri4 <- A[, 21]

### colors of weather and yield time series lines

library(RColorBrewer)
library(ggplot2)
display.brewer.all()

## yield 
barplot(rep(1,9),col=brewer.pal(9,"GnBu")[1:9])
myColors_yield=brewer.pal(9,"GnBu")[1:9] 
barplot(rep(1,9),col= myColors_yield, main="GnBu")
n=96
barplot(rep(1, n),col= colorRampPalette(colors = myColors_yield)(n),main="spline")
colorRampPalette(colors = myColors_yield)(n)

## detrended yield
barplot(rep(1,9),col=brewer.pal(9,"Purples")[1:9])
myColors_deyield=brewer.pal(9,"Purples")[1:9] 
barplot(rep(1,9),col= myColors_deyield, main="Purples")
n=96
barplot(rep(1, n),col= colorRampPalette(colors = myColors_deyield)(n),main="spline")
colorRampPalette(colors = myColors_deyield)(n)

## yield loss 
barplot(rep(1,9),col=brewer.pal(9,"RdPu")[1:9])
myColors_yieldlo=brewer.pal(9,"RdPu")[1:9] 
barplot(rep(1,9),col= myColors_yieldlo, main="RdPu")
n=96
barplot(rep(1, n),col= colorRampPalette(colors = myColors_yieldlo)(n),main="spline")
colorRampPalette(colors = myColors_yieldlo)(n)

## gdd
barplot(rep(1,9),col=brewer.pal(9,"OrRd")[1:9])
myColors_gdd=brewer.pal(9,"OrRd")[1:9] 
barplot(rep(1,9),col= myColors_gdd, main="OrRd")
n=96
barplot(rep(1, n),col= colorRampPalette(colors = myColors_gdd)(n),main="spline")
colorRampPalette(colors = myColors_gdd)(n)

## rainfall
barplot(rep(1,9),col=brewer.pal(9,"Greens")[1:9])
myColors_rainfall=brewer.pal(9,"Greens")[1:9] 
barplot(rep(1,9),col= myColors_rainfall, main="Greens")
n=96
barplot(rep(1, n),col= colorRampPalette(colors = myColors_rainfall)(n),main="spline")
colorRampPalette(colors = myColors_rainfall)(n)


### plot yield, detrended yield and yield loss
range(yield)
range(detyield)
range(yieldloss)

library(grid)
library(gridExtra)
# define a file path and a file name here.
pdf(file="C\\...\\filename.pdf", height=18, width=14)
par(mfrow = c(3, 1))

plot(year[1:23],yield[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(15, 81), col=alpha(colorRampPalette(colors = myColors_yield)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b],yield[a:b], col=alpha(colorRampPalette(colors = myColors_yield)(n)[i],0.8))
} 
title(xlab="Year", cex.lab = 2)
title(ylab = "Yield (bushels/acre)", cex.lab = 2, line = 2.5)

plot(year[1:23],detyield[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(16, 57), col=alpha(colorRampPalette(colors = myColors_deyield)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b],detyield[a:b], col=alpha(colorRampPalette(colors = myColors_deyield)(n)[i],0.8))
}
title(xlab="Year", cex.lab = 2)
title(ylab = "Detrended Yield (bushels/acre)", cex.lab = 2, line = 2.5)

plot(year[1:23],yieldloss[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(0, 36), col=alpha(colorRampPalette(colors = myColors_yieldlo)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b],yieldloss[a:b], col=alpha(colorRampPalette(colors = myColors_yieldlo)(n)[i],0.8))
}
title(xlab="Year", cex.lab = 2)
title(ylab = "Yield Loss (bushels/acre)", cex.lab = 2, line = 2.5)

dev.off()


####### plot whole-cycle GDD and CRI
range(gdd0)
range(cri0)

library(grid)
library(gridExtra)
# define a file path and file name here.
pdf(file="C:\\...\\filename.pdf", height=8, width=20)
par(mfrow = c(1, 2))
plot(year[1:23],gdd0[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(1100, 2400), col=alpha(colorRampPalette(colors = myColors_gdd)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b],gdd0[a:b], col=alpha(colorRampPalette(colors = myColors_gdd)(n)[i],0.8))
}
title(xlab="Year", cex.lab = 2)
title(ylab = "GDD", cex.lab = 2, line = 2.5)

plot(year[1:23],cri0[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(220, 1050), col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b],cri0[a:b], col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[i],0.8))
}
title(xlab="Year", cex.lab = 2)
title(ylab = "CR", cex.lab = 2, line = 2.5)
dev.off()

## plot phase-division GDD 
range(gdd1)
range(gdd2)
range(gdd3)
range(gdd4)

library(grid)
library(gridExtra)
# define a fine path and a file name here.
pdf(file="C:\\...\\filename.pdf", height=8, width=10)
par(mfrow = c(2, 2))

plot(year[1:23],gdd1[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(150, 1050), col=alpha(colorRampPalette(colors = myColors_gdd)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b],gdd1[a:b], col=alpha(colorRampPalette(colors = myColors_gdd)(n)[i],0.8))
}
title(main="Emerged", cex.main = 2, adj = 0)
title(xlab="Year", cex.lab = 2)
title(ylab = "GDD", cex.lab = 2, line = 2.5)

plot(year[1:23],gdd2[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(150, 1050), col=alpha(colorRampPalette(colors = myColors_gdd)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b],gdd2[a:b], col=alpha(colorRampPalette(colors = myColors_gdd)(n)[i],0.8))
}
title(main="Blooming", cex.main = 2, adj = 0)
title(xlab="Year", cex.lab = 2)
title(ylab = "GDD", cex.lab = 2, line = 2.5)


plot(year[1:23],gdd3[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(150, 1050), col=alpha(colorRampPalette(colors = myColors_gdd)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b],gdd3[a:b], col=alpha(colorRampPalette(colors = myColors_gdd)(n)[i],0.8))
}
title(main="Setting pods", cex.main = 2, adj = 0)
title(xlab="Year", cex.lab = 2)
title(ylab = "GDD", cex.lab = 2, line = 2.5)

plot(year[1:23],gdd4[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(150, 1050), col=alpha(colorRampPalette(colors = myColors_gdd)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b],gdd4[a:b], col=alpha(colorRampPalette(colors = myColors_gdd)(n)[i],0.8))
}
title(main="Dropping leaves", cex.main = 2, adj = 0)
title(xlab="Year", cex.lab = 2)
title(ylab = "GDD", cex.lab = 2, line = 2.5)

dev.off()

## plot phase-division CRI

range(cri1)
range(cri2)
range(cri3)
range(cri4)

library(grid)
library(gridExtra)
# define a fine path and a file name here.
pdf(file="C:\\...\\filename.pdf", height=8, width=10)
par(mfrow = c(2, 2))

plot(year[1:23], cri1[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(15, 550), col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b], cri1[a:b], col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[i],0.8))
}
title(main="Emerged", cex.main = 2, adj = 0)
title(xlab="Year", cex.lab = 2)
title(ylab = "CR", cex.lab = 2, line = 2.5)

plot(year[1:23], cri2[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(15, 550), col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b], cri2[a:b], col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[i],0.8))
}
title(main="Blooming", cex.main = 2, adj = 0)
title(xlab="Year", cex.lab = 2)
title(ylab = "CR", cex.lab = 2, line = 2.5)


plot(year[1:23], cri3[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(15, 550), col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b], cri3[a:b], col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[i],0.8))
}
title(main="Setting pods", cex.main = 2, adj = 0)
title(xlab="Year", cex.lab = 2)
title(ylab = "CR", cex.lab = 2, line = 2.5)

plot(year[1:23], cri4[1:23],xlab="", ylab="", cex.axis = 1.8, type="l", ylim=c(15, 550), col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[1],0.8))
for (i in 2:96){
  a = (i-1)*23+1
  b = 23*i
  lines(year[a:b], cri4[a:b], col=alpha(colorRampPalette(colors = myColors_rainfall)(n)[i],0.8))
}
title(main="Dropping leaves", cex.main = 2, adj = 0)
title(xlab="Year", cex.lab = 2)
title(ylab = "CR", cex.lab = 2, line = 2.5)

dev.off()


