rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rio)

dataT <- import("Train16.xlsx")

A <- as.matrix(dataT[, 4:25])
class(A) <- "numeric"

gdd0T <- A[, 8]
cri0T <- A[, 9]
#rdi0T <- A[, 10]

# quadratic utility
yieldlossT <- A[, 7]

# data scaling
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
yT <- range01(yieldlossT)
x01T <- range01(gdd0T)
x02T <-range01(cri0T)
#x02T <- range01(rdi0T)

# exponential utility
#ALPHA <- c(0.0052, 0.008, 0.0103)
#alpha <- ALPHA[1]
#ylTexp <- exp(alpha*yieldlossT)
#yeT <- range01(ylTexp)

## Validation dataset

dataV <- import("Validation4.xlsx")

B <- as.matrix(dataV[, 4:25])
class(B) <- "numeric"

gdd0V <- B[, 8]
cri0V <- B[, 9]

### quadratic utility
yieldlossV <- B[, 7]

# data scaling
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
yV <- range01(yieldlossV)
x01V <- range01(gdd0V)
x02V <-range01(cri0V)

# exponential utility
#ALPHA <- c(0.0052, 0.008, 0.0103)
#alpha <- ALPHA[1]
#ylVexp <- exp(alpha*yieldlossV)
#yeV <- range01(ylVexp)

source("basic-functions.R")

source("Train-W.R")

source("Validation-W.R")

nsegseq <- seq(21, 39, by = 3)
length(nsegseq)

# offset vectors and sequence
M <- PSA <- b.seq <- PRED <- list()
it.seq <- dla.seq <- c.time <- mse.seq <- mse.new.seq <- rmse.seq <- rmse.new.seq <- ratio.seq <- vector(l = length(nsegseq))

for (i in (1:length(nsegseq))){
  start <- proc.time()[3]
  M[[i]] <- NEWXZ(x01V, x02V, nsegseq[i])
  PSA[[i]] <- WCPSA(yT, x01T, x02T, nsegseq[i])
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

