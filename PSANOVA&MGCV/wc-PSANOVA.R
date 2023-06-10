##############################################
# File: wc-PSANOVA.r 
# Author: [Jing Zou] 
# Created: [2023.06.10] 
# Description: [estimate whole-cycle models using PSANOVA] 
##############################################

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

source("basic-functions.R")

source("Train-W.R")

###### quadratic utility case

WC <- WCPSA(y, x01, x02, 21)

Fit_wc_back <- back(WC$Fit, yieldloss)

# compare the range of fitted yieldloss and yieldloss
range(Fit_wc_back)
range(yieldloss)

# RMSE and adj.R^2
y_back <- back(y, yieldloss)
RMSE_wc <- sqrt(mean((y_back - Fit_wc_back)^2))
RMSE_wc

expvar1 <- sum((y-WC$Fit)^2)
expvar2 <- sum((y-mean(y))^2)
R2wc = 1 - expvar1/expvar2
n = length(y)
k = 2
adR2wc = 1- (1-R2wc)*(n-1)/(n-k-1)
adR2wc

###### exponential utility case 
## under exponential utility, reponse is ylexp
ALPHA <- c(0.0052, 0.008, 0.0103)
## take the high level risk aversion for example
alpha <- ALPHA[3]
ylexp <- exp(alpha*yieldloss)
ye <- range01(ylexp)

WCe <- WCPSA(ye, x01, x02, 30)

## RMSE
ye_back <- back(ye,ylexp)
y_exp <- (1/alpha)*log(ye_back)
Fit_wce_back <- back(WCe$Fit, ylexp)
Fit_wce <- (1/alpha)*log(Fit_wce_back)
RMSE_wce <- sqrt(mean((y_exp - Fit_wce)^2))
RMSE_wce

## adj.R^2 is scale-independent.
expvar1e <- sum((ye-WCe$Fit)^2)
expvar2 <- sum((ye-mean(ye))^2)
R2wce = 1 - expvar1e/expvar2
n = length(ye)
k = 2
adR2wce = 1- (1-R2wce)*(n-1)/(n-k-1)
adR2wce

###############Visualization with contour plots################# 
library(ggplot2)
library(akima)
library(reshape2)
library(metR)
library(grid)
library(gridExtra)

xx1 <- gdd0
xx2 <- cri0
z1 = Fit_wc_back

yp1 <- as.numeric(z1)
dfp1 <- data.frame(cbind(xx1,xx2,yp1))

dfp1.int <- interp(x = dfp1$xx1, y = dfp1$xx2, z=dfp1$yp1,
                   xo=seq(min(dfp1$xx1), max(dfp1$xx1), length = 2208),
                   yo=seq(min(dfp1$xx2), max(dfp1$xx2), length = 2208))

datp1 <- melt(dfp1.int$z, na.rm = TRUE)
names(datp1) <- c("GDD", "Rainfall", "f1")
datp1$GDD <- dfp1.int$x[datp1$GDD]
datp1$Rainfall <- dfp1.int$y[datp1$Rainfall]

#library(graphics)
#library(exams)
range(gdd0)
range(cri0)
range(z1)

p <- ggplot(datp1, aes(GDD, Rainfall, z = f1)) + 
  geom_contour(aes(z = f1, colour=stat(level)), linewidth=2) + 
  geom_text_contour(aes(z = f1), stroke = 0.5, label.placer = label_placer_flattest(), check_overlap = TRUE, size = 7) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = "transparent", colour = "transparent")) + 
  scale_x_continuous(limits = c(1000, 2500), n.breaks = 10) +
  scale_y_continuous(limits = c(200, 1100), n.breaks = 10) +
  labs(colour = "Yield Loss") +
  scale_color_viridis_c(option="B") +
  expand_limits(colour = seq(0,30.5)) +
  labs(x=expression('GDD'[0]), y=expression('CR'[0]))+
  theme(text = element_text(size=24))

## define a file path and filename.pdf
pdf(file="C:\\...\\filename.pdf", height=16, width=24)
p
dev.off()
