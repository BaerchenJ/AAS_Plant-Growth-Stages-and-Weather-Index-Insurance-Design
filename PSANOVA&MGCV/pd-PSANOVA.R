
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


######################### PSANOVA ###################

source("basic-functions.R")

source("Train-P.R")

###### quadratic utility case

PD <- PDPSA(y, x11, x12, x21, x22, x31, x32, x41, x42, 39)

Fit_pd_back <- back(PD$Fit, yieldloss)

# compare the range of fitted yieldloss and yieldloss
range(Fit_pd_back)
range(yieldloss)

# RMSE and adj.R^2
y_back <- back(y, yieldloss)
RMSE_pd <- sqrt(mean((y_back - Fit_pd_back)^2))
RMSE_pd

expvar1 <- sum((y-PD$Fit)^2)
expvar2 <- sum((y-mean(y))^2)
R2pd = 1 - expvar1/expvar2
n = length(y)
k = 8
adR2pd = 1- (1-R2pd)*(n-1)/(n-k-1)
adR2pd

###### exponential utility case 
## under exponential utility, reponse is ylexp
ALPHA <- c(0.0052, 0.008, 0.0103)
## take the low level risk aversion for example
alpha <- ALPHA[1]
ylexp <- exp(alpha*yieldloss)
ye <- range01(ylexp)

PDe <- PDPSA(ye, x11, x12, x21, x22, x31, x32, x41, x42, 39)

## RMSE
ye_back <- back(ye,ylexp)
y_exp <- (1/alpha)*log(ye_back)
Fit_pde_back <- back(PDe$Fit, ylexp)
Fit_pde <- (1/alpha)*log(Fit_pde_back)
RMSE_pde <- sqrt(mean((y_exp - Fit_pde)^2))
RMSE_pde

## adj.R^2 is scale-independent.
expvar1e <- sum((ye-PDe$Fit)^2)
expvar2 <- sum((ye-mean(ye))^2)
R2pde = 1 - expvar1e/expvar2
n = length(ye)
k = 8
adR2pde = 1- (1-R2pde)*(n-1)/(n-k-1)
adR2pde

############################### Visualization ###############################
####################### Yield loss contribution extraction ##################

idx <- PD$idx
M <- PD$M
b <- PD$b

# Stage 1
index1 <- 1 * (idx == 5)
index1[1] <- 1
index1[2] <- 1
Fit1 <- as.numeric(M %*% (index1 * b))

index2 <- 1 * (idx == 6)
index2[3] <- 1
Fit2 <- as.numeric(M %*% (index2 * b))

index12a <- 1 * (idx == 7)
Fit12a <- as.numeric(M %*% (index12a * b))

index12b <- 1 * (idx == 8)
Fit12b <- as.numeric(M %*% (index12b * b))

index12c <- 1 * (idx == 9)
index12c[4] <- 1
Fit12c <- as.numeric(M %*% (index12c * b))

FitS1 <- Fit1 + Fit2 + Fit12a + Fit12b + Fit12c
FitS1_back <- back(FitS1, yieldloss)
Yield_Loss_1 <- FitS1_back

# Stage 2
index3 <- 1 * (idx == 10)
index3[5] <- 1
index3[6] <- 1
Fit3 <- as.numeric(M %*% (index3 * b))

index4 <- 1 * (idx == 11)
index4[7] <- 1
Fit4 <- as.numeric(M %*% (index4 * b))

index34a <- 1 * (idx == 12)
Fit34a <- as.numeric(M %*% (index34a * b))

index34b <- 1 * (idx == 13)
Fit34b <- as.numeric(M %*% (index34b * b))

index34c <- 1 * (idx == 14)
index34c[8] <- 1
Fit34c <- as.numeric(M %*% (index34c * b))

FitS2 <- Fit3 + Fit4 + Fit34a + Fit34b + Fit34c
FitS2_back <- back(FitS2, yieldloss)
Yield_Loss_2 <- FitS2_back

# Stage 3
index5 <- 1 * (idx == 15)
index5[9] <- 1
index5[10] <- 1
Fit5 <- as.numeric(M %*% (index5 * b))

index6 <- 1 * (idx == 16)
index6[11] <- 1
Fit6 <- as.numeric(M %*% (index6 * b))

index56a <- 1 * (idx == 17)
Fit56a <- as.numeric(M %*% (index56a * b))

index56b <- 1 * (idx == 18)
Fit56b <- as.numeric(M %*% (index56b * b))

index56c <- 1 * (idx == 19)
index56c[12] <- 1
Fit56c <- as.numeric(M %*% (index56c * b))

FitS3 <- Fit5 + Fit6 + Fit56a + Fit56b + Fit56c
FitS3_back <- back(FitS3, yieldloss)
Yield_Loss_3 <- FitS3_back

## Stage 4
index7 <- 1 * (idx == 20)
index7[13] <- 1
index7[14] <- 1
Fit7 <- as.numeric(M %*% (index7 * b))

index8 <- 1 * (idx == 21)
index8[15] <- 1
Fit8 <- as.numeric(M %*% (index8 * b))

index78a <- 1 * (idx == 22)
Fit78a <- as.numeric(M %*% (index78a * b))

index78b <- 1 * (idx == 23)
Fit78b <- as.numeric(M %*% (index78b * b))

index78c <- 1 * (idx == 24)
index78c[16] <- 1
Fit78c <- as.numeric(M %*% (index78c * b))

FitS4 <- Fit7 + Fit8 + Fit78a + Fit78b + Fit78c
FitS4_back <- back(FitS4, yieldloss)
Yield_Loss_4 <- FitS4_back


z1 <- Yield_Loss_1
z2 <- Yield_Loss_2
z3 <- Yield_Loss_3
z4 <- Yield_Loss_4


###### Visualization ####### contour plots ########

library(ggplot2)
library(akima)
library(reshape2)
library(metR)

## data transformation for contour plots

# Stage 1
xx1 <- gdd1
xx2 <- cri1

yp1 <- as.numeric(z1)
dfp1 <- data.frame(cbind(xx1,xx2,yp1))

dfp1.int <- interp(x = dfp1$xx1, y = dfp1$xx2, z=dfp1$yp1,
                   xo=seq(min(dfp1$xx1), max(dfp1$xx1), length = 2208),
                   yo=seq(min(dfp1$xx2), max(dfp1$xx2), length = 2208))

datp1 <- melt(dfp1.int$z, na.rm = TRUE)
names(datp1) <- c("GDD1", "Rainfall1", "f.1")
datp1$GDD1 <- dfp1.int$x[datp1$GDD1]
datp1$Rainfall1 <- dfp1.int$y[datp1$Rainfall1]

# Stage 2
xx3 <- gdd2
xx4 <- cri2

yp2 <- as.numeric(z2)
dfp2 <- data.frame(cbind(xx3,xx4,yp2))
dfp2.int <- interp(x = dfp2$xx3, y = dfp2$xx4, z=dfp2$yp2,
                   xo=seq(min(dfp2$xx3), max(dfp2$xx3), length = 2208),
                   yo=seq(min(dfp2$xx4), max(dfp2$xx4), length = 2208))

datp2 <- melt(dfp2.int$z, na.rm = TRUE)
names(datp2) <- c("GDD2", "Rainfall2", "f.2")
datp2$GDD2 <- dfp2.int$x[datp2$GDD2]
datp2$Rainfall2 <- dfp2.int$y[datp2$Rainfall2]

# Stage 3
xx5 <- gdd3
xx6 <- cri3

yp3 <- as.numeric(z3)
dfp3 <- data.frame(cbind(xx5,xx6,yp3))
dfp3.int <- interp(x = dfp3$xx5, y = dfp3$xx6, z=dfp3$yp3,
                   xo=seq(min(dfp3$xx5), max(dfp3$xx5), length = 2208),
                   yo=seq(min(dfp3$xx6), max(dfp3$xx6), length = 2208))

datp3 <- melt(dfp3.int$z, na.rm = TRUE)
names(datp3) <- c("GDD3", "Rainfall3", "f.3")
datp3$GDD3 <- dfp3.int$x[datp3$GDD3]
datp3$Rainfall3 <- dfp3.int$y[datp3$Rainfall3]

# Stage 4
xx7 <- gdd4
xx8 <- cri4

yp4 <- as.numeric(z4)
dfp4 <- data.frame(cbind(xx7,xx8,yp4))
dfp4.int <- interp(x = dfp4$xx7, y = dfp4$xx8, z=dfp4$yp4,
                   xo=seq(min(dfp4$xx7), max(dfp4$xx7), length = 2208),
                   yo=seq(min(dfp4$xx8), max(dfp4$xx8), length = 2208))

datp4 <- melt(dfp4.int$z, na.rm = TRUE)
names(datp4) <- c("GDD4", "Rainfall4", "f.4")
datp4$GDD4 <- dfp4.int$x[datp4$GDD4]
datp4$Rainfall4 <- dfp4.int$y[datp4$Rainfall4]

## define the domain shown by xaxis and yaxis in plots
## based on the ranges of yield loss contributions and weather indices
# range of yield loss contribution
min(cbind(z1, z2, z3, z4))
max(cbind(z1, z2, z3, z4))

# range of GDD
min(cbind(gdd1, gdd2, gdd3, gdd4))
max(cbind(gdd1, gdd2, gdd3, gdd4))

# range of CRI
min(cbind(cri1, cri2, cri3, cri4))
max(cbind(cri1, cri2, cri3, cri4))


## Plotting now the phase-division contour plots

p1 <- ggplot(datp1, aes(GDD1, Rainfall1, z = f.1)) + 
  geom_contour(aes(z = f.1, colour=stat(level)), size=2) + 
  geom_text_contour(aes(z = f.1), stroke = 0.5, label.placer = label_placer_flattest(), check_overlap = TRUE, size = 5) + 
  ggtitle("Emerged") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = "transparent", colour = "transparent")) + 
  scale_x_continuous(limits = c(150, 1050), n.breaks = 10) +
  scale_y_continuous(limits = c(15, 535), n.breaks = 10) +
  scale_color_viridis_c() +
  expand_limits(colour = seq(-10, 16)) +
  labs(x=expression('GDD'[1]), y=expression('CR'[1]))+
  theme(text = element_text(size=20)) +
  labs(colour = "Yield Loss")

p2 <- ggplot(datp2, aes(GDD2, Rainfall2, z = f.2)) + 
  geom_contour(aes(z = f.2, colour=stat(level)), size=2) + 
  geom_text_contour(aes(z = f.2), stroke = 0.5, label.placer = label_placer_flattest(), check_overlap = TRUE, size = 5) + 
  ggtitle("Blooming") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = "transparent", colour = "transparent")) + 
  scale_x_continuous(limits = c(150, 1050), n.breaks = 10) +
  scale_y_continuous(limits = c(15, 535), n.breaks = 10) + 
  scale_color_viridis_c() +
  expand_limits(colour = seq(-10, 16)) +
  labs(x=expression('GDD'[2]), y=expression('CR'[2]))+
  theme(text = element_text(size=20)) +
  labs(colour = "Yield Loss")

p3 <- ggplot(datp3, aes(GDD3, Rainfall3, z = f.3)) + 
  geom_contour(aes(z = f.3, colour=stat(level)), size=2) + 
  geom_text_contour(aes(z = f.3), stroke = 0.5, label.placer = label_placer_flattest(), check_overlap = TRUE, size = 5) + 
  ggtitle("Setting pods") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill = "transparent", colour = "transparent")) + 
  scale_x_continuous(limits = c(150, 1050), n.breaks = 10) +
  scale_y_continuous(limits = c(15, 535), n.breaks = 10) +
  scale_color_viridis_c() +
  expand_limits(colour = seq(-10, 16)) +
  labs(x=expression('GDD'[3]), y=expression('CR'[3])) +
  theme(text = element_text(size=20)) +
  labs(colour = "Yield Loss")

p4 <- ggplot(datp4, aes(GDD4, Rainfall4, z = f.4)) + 
  geom_contour(aes(z = f.4, colour=stat(level)), size=2) + 
  geom_text_contour(aes(z = f.4), stroke = 0.5, label.placer = label_placer_flattest(), check_overlap = TRUE, size = 5) + 
  ggtitle("Dropping leaves") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.key = element_rect(fill = "transparent", colour = "transparent")) + 
  scale_x_continuous(limits = c(150, 1050), n.breaks = 10) +
  scale_y_continuous(limits = c(15, 535), n.breaks = 10) +
  scale_color_viridis_c() +
  expand_limits(colour = seq(-10, 16)) +
  labs(x=expression('GDD'[4]), y=expression('CR'[4])) +
  theme(text = element_text(size=20)) +
  labs(colour = "Yield Loss")

## save phase-division plot in pdf. file

library(grid)
library(gridExtra)
## define the file path and a "filename" in pdf() please
pdf(file="C:\\...\\filename.pdf", height=18, width=26)
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol=2)
dev.off()

