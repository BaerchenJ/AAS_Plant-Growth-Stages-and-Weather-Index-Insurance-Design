rm(list=ls())

library(ncdf4)
library(pracma)

## define a file path and a file name
prcp20 <- nc_open("D://...//prcp_2020_illi.nc")

## define the list of Rainfall data in each Stage 
prS1 <- prS2 <- prS3 <- prS4 <- list()

for (i in (124:180)){
  prS1[[length(prS1)+1]] <- ncvar_get(prcp20,varid=strcat("Band", num2str(i,'%d')))
}

for (i in (173:229)){
  prS2[[length(prS2)+1]] <- ncvar_get(prcp20,varid=strcat("Band", num2str(i,'%d')))
}

for (i in (187:243)){
  prS3[[length(prS3)+1]] <- ncvar_get(prcp20,varid=strcat("Band", num2str(i,'%d')))
}

for (i in (250:292)){
  prS4[[length(prS4)+1]] <- ncvar_get(prcp20,varid=strcat("Band", num2str(i,'%d')))
}

## weekly precip in each phase

prS1WA <- list()

aa = (length(prS1)-1)/7

for (i in 1 : aa){
  bb = (i-1)*7+1
  cc = i*7
  prS1WA[[i]] <- (1/7) * Reduce("+", prS1[c(bb : cc)])
}

prS2WA <- list()

aa = (length(prS2)-1)/7

for (i in 1 : aa){
  bb = (i-1)*7+1
  cc = i*7
  prS2WA[[i]] <- (1/7) * Reduce("+", prS2[c(bb : cc)])
}

prS3WA <- list()

aa = (length(prS3)-1)/7

for (i in 1 : aa){
  bb = (i-1)*7+1
  cc = i*7
  prS3WA[[i]] <- (1/7) * Reduce("+", prS3[c(bb : cc)])
}

prS4WA <- list()

aa = (length(prS4)-1)/7

for (i in 1 : aa){
  bb = (i-1)*7+1
  cc = i*7
  prS4WA[[i]] <- (1/7) * Reduce("+", prS4[c(bb : cc)])
}

#### Rainfall Deficit Index

RDI20S3_list <- list()

Daily20S3_list <- list()

for (j in 1 : 8){
aa = (j-1)*7+1
bb = j*7

for (i in aa : bb){
Daily20S3_list[[i]] <- prS3[[i]] - prS3WA[[j]]
}

M <- Reduce("+", Daily20S3_list)
MM <- matrix(0, nrow = 347, ncol = 596)

RDI20S3_list[[j]] <- pmin(MM, M)
}


RDI20S3 <- Reduce("+", RDI20S3_list)


lon <- ncdim_def("longitude", "degrees", prcp20$dim$x$vals)
lat <- ncdim_def("latitude", "degrees", prcp20$dim$y$vals)

mv <- -9999
rdi20s <- ncvar_def("precipitation", "mm", list(lon,lat), mv)

setwd("D:\\EUO\\Dataset\\Illinois\\IllinoisSOYBEAN\\RDI")

ncnew <- nc_create("rdi20s.nc", list(rdi20s))
ncvar_put(ncnew, rdi20s, RDI20S3)
nc_close(ncnew)

## after obtaining this new nc.file with RDI, 
## we extract the average values in ArcGIS
## according to county polygons.


## Finally, we link all the RDI with county and year
## in this file "WeatherIndices&Yields.xlsx"
