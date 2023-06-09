
rm(list=ls())

source("basic-functions.r")

source("Train-P.r")

source("Validation-P.r")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rio)

dataT <- import("Train16.xlsx")

A <- as.matrix(dataT[, 4:25])
class(A) <- "numeric"

## Stage 1
gdd1T <- A[, 11]
cri1T <- A[, 12]
#rdi1T <- A[, 13]

# Stage 2
gdd2T <- A[, 14]
cri2T <- A[, 15]
#rdi2T <- A[, 16]

# Stage 3 
gdd3T <- A[, 17]
cri3T <- A[, 18]
#rdi3T <- A[, 19]

# Stage 4 
gdd4T <- A[, 20]
cri4T <- A[, 21]
#rdi4T <- A[, 22]

yieldlossT <- A[, 7]

# data scaling
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
yT <- range01(yieldlossT)

# phase-division
x11T <- range01(gdd1T)
x12T <- range01(cri1T)
# x12T <- range01(rdi1T)

x21T <- range01(gdd2T)
x22T <- range01(cri2T)
# x22T <- range01(rdi2T)

x31T <- range01(gdd3T)
x32T <- range01(cri3T)
# x32T <- range01(rdi3T)

x41T <- range01(gdd4T)
x42T <- range01(cri4T)
# x42T <- range01(rdi4T)

# exponential utility
#ALPHA <- c(0.0052, 0.008, 0.0103)
#alpha <- ALPHA[1]
#ylTexp <- exp(alpha*yieldlossT)
#yeT <- range01(ylTexp)

## Validation dataset

dataV <- import("Validation4.xlsx")

B <- as.matrix(dataV[, 4:25])
class(B) <- "numeric"

## Stage 1
gdd1V <- B[, 11]
cri1V <- B[, 12]
#rdi1V <- B[, 13]

# Stage 2
gdd2V <- B[, 14]
cri2V <- B[, 15]
#rdi2V <- B[, 16]

# Stage 3 
gdd3V <- B[, 17]
cri3V <- B[, 18]
#rdi3V <- B[, 19]

# Stage 4 
gdd4V <- B[, 20]
cri4V <- B[, 21]
#rdi4V <- B[, 22]

yieldlossV <- B[, 7]

# data scaling
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
yV <- range01(yieldlossV)

# phase-division
x11V <- range01(gdd1V)
x12V <- range01(cri1V)
# x12V <- range01(rdi1V)

x21V <- range01(gdd2V)
x22V <- range01(cri2V)
# x22V <- range01(rdi2V)

x31V <- range01(gdd3V)
x32V <- range01(cri3V)
# x32V <- range01(rdi3V)

x41V <- range01(gdd4V)
x42V <- range01(cri4V)
# x42V <- range01(rdi4V)

# exponential utility
#ALPHA <- c(0.0052, 0.008, 0.0103)
#alpha <- ALPHA[1]
#ylVexp <- exp(alpha*yieldlossV)
#yeV <- range01(ylVexp)

nsegseq <- seq(21, 39, by = 3)
length(nsegseq)

# offset vectors and sequence
M <- PSA <- b.seq <- PRED <- list()
it.seq <- dla.seq <- c.time <- mse.seq <- mse.new.seq <- rmse.seq <- rmse.new.seq <- ratio.seq <- vector(l = length(nsegseq))

for (i in (1:length(nsegseq))){
start <- proc.time()[3]
M[[i]] <- NEWXZ(x11V, x12V, x21V, x22V, x31V, x32V, x41V, x42V, nsegseq[i])
PSA[[i]] <- PDPSA(yT, x11T, x12T, x21T, x22T, x31T, x32T, x41T, x42T, nsegseq[i])
end <- proc.time()[3]
it.seq[i] <- PSA[[i]]$it
dla.seq[i] <- PSA[[i]]$dla
c.time[i] <- end - start
b.seq[[i]] <- PSA[[i]]$b
mse.seq[i] <- PSA[[i]]$mse
PRED[[i]] <- M[[i]]%*%b.seq[[i]]
mse.new.seq[i] <- mean((yV-PRED[[i]])^2)
rmse.seq[i] <- sqrt(mse.seq[i])
rmse.new.seq[i] <-sqrt(mse.new.seq[i])
ratio.seq[i] <- rmse.new.seq[i]/rmse.seq[i]
res <- list(it.seq = it.seq, dla.seq = dla.seq, ratio.seq = ratio.seq, rmse.new.seq = rmse.new.seq, rmse.seq = rmse.seq, c.time = c.time)
res
}

distance <- dla.seq - 1e-06
sel.seq <- cbind(nsegseq, it.seq, distance, ratio.seq, rmse.seq, rmse.new.seq, c.time)

library(writexl)
sn <- as.data.frame(sel.seq)
## define a file path and a file name here to save the results.
write_xlsx(sn, "C:\\...\\filename.xlsx")
