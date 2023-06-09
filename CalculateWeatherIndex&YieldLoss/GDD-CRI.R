rm(list=ls())

library(ncdf4)
## firstly download the nc.file from Daymet dataset 
## and define a file path
## e.g., we download the 2021 precipitation nc.file 
## then, we cut it in QGIS according to the shape file of Illinois State,
## finally, we obtain the file "prcp_2021_illi.nc".

prcp21 <- nc_open("D://...//prcp_2021_illi.nc")

A <- ncvar_get(prcp21,"Band1")
str(A)

library(pracma)
M <- A-A 

## the range of i is determined by start DOY and end DOY
## e.g., the start and end date of the Stage 4 in 2021
## for Illinois Soybeans are 241 and 290.

for (i in (241:290)){
  M <- ncvar_get(prcp21,varid=strcat("Band", num2str(i,'%d'))) + M
}
print(M)

lon <- ncdim_def("longitude", "degrees", prcp21$dim$x$vals)
lat <- ncdim_def("latitude", "degrees", prcp21$dim$y$vals)

mv <- -9999
prcp21d <- ncvar_def("precipitation", "mm", list(lon,lat), mv)

ncnew <- nc_create("prcp21d.nc", list(prcp21d))

ncvar_put(ncnew, prcp21d, M)
nc_close(ncnew)

## after obtaining this new nc.file with CRI, 
## we extract the average values in ArcGIS
## according to county polygons.

#########################################################################################################

## GDD Calculation
rm(list=ls())

library(ncdf4)
library(pracma)

## nc.file downloading and pre-processing is similar with precipitation ones.
tmax20 <- nc_open("D://...//tmax_2020_illi.nc")
tmin20 <- nc_open("D://...//tmin_2020_illi.nc")

## sdat, edat are start and end dates of each phase
## tbase is the baseline temperature, i.e., 10 degree Celcius in our study.
SGDD <- function(sdat, edat, tbase){
    nc = edat-sdat+1
    TMAX <- matrix(NA, 347, nc*596)
    TMIN <- matrix(NA, 347, nc*596)
    TAVE <- matrix(NA, 347, nc*596)
    SGDD <- matrix(NA, 347, nc*596)
    for (i in (sdat:edat)){
    a = (i-sdat)*596+1
    b = (i-sdat+1)*596
    TMAX[,a:b] <- ncvar_get(tmax20, varid=strcat("Band", num2str(i,'%d')))
    TMIN[,a:b] <- ncvar_get(tmin20, varid=strcat("Band", num2str(i,'%d')))
    TAVE[,a:b] = (TMAX[,a:b] + TMIN[,a:b])/2
    SGDD[,a:b] <- pmax((TAVE[,a:b] - tbase), 0)
    }
    return(SGDD)
}

A <- SGDD(124,292,10)

###################### create a list 
mylist <- list()
c = sdat = 124
d = edat = 292
for (i in (c:d)){
  a = (i-c)*596+1
  b = (i-c+1)*596
  mylist[[i-c+1]] <- A[,a:b] 
}

###################### sum as GDD
M <- Reduce("+",mylist)

lon <- ncdim_def("longitude", "degrees", tmax20$dim$x$vals)
lat <- ncdim_def("latitude", "degrees", tmax20$dim$y$vals)

mv <- -9999
gdd20 <- ncvar_def("GDD", "celcius", list(lon,lat), mv)

ncnew <- nc_create("gdd20.nc", list(gdd20))

ncvar_put(ncnew, gdd20, M)
nc_close(ncnew)

## after obtaining this new nc.file with GDD, 
## we extract the average values in ArcGIS
## according to county polygons.


## Finally, we link all the GDD and CRI with county and year
## in this file "WeatherIndices&Yields.xlsx"


