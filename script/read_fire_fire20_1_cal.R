rm(list = ls())
graphics.off()
gc()

require(rgdal)
library(R.matlab)
library(raster)
library(rgeos)
library(secr)


dir_fire<-"/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/"

# The input file geodatabase
mtbs <- readOGR("/Users/marco/Documents/dati/fire_us/fire20_1_california.shp", verbose = FALSE)
summary(mtbs)
mtbs = spTransform(mtbs,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))




#plot(rr)

BA = array(data = 0, dim = c(1, length(1950:2020) * 12))

anni=as.numeric(format(as.Date(mtbs@data$ALARM_D),"%Y"))
mesi=as.numeric(format(as.Date(mtbs@data$ALARM_D),"%m"))

k = 0
for (iyear in 1950:2020) {
  for (imonth in 1:12) {
    k = k + 1
    idx2 = which(anni == iyear &
                   mesi == imonth)
    
    
      if (length(idx2) >= 1) {
        BA[1, k] = sum(mtbs[idx2,]@data$GIS_ACR)
      }
      
  }
}



save(BA, file = file.path(dir_fire,"fire_1950_2020_monthly_cal.RData"))
writeMat(file.path(dir_fire, paste("fire_1950_2020_monthly_cal.mat",sep="")), BA=BA)

