rm(list = ls())
graphics.off()
gc()

require(rgdal)
library(R.matlab)
library(raster)
library(rgeos)
library(secr)
library(pracma)

dir_fire<-"/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/"

# The input file geodatabase
mtbs <- readOGR("/Users/marco/Documents/dati/fire_us/fire20_1_baileys.shp", verbose = FALSE)
summary(mtbs)
mtbs = spTransform(mtbs,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))




#plot(rr)

BA = array(data = 0, dim = c(1, length(1950:2020) * 12))

anni=as.numeric(format(as.Date(mtbs@data$ALARM_D),"%Y"))
mesi=as.numeric(format(as.Date(mtbs@data$ALARM_D),"%m"))
giorni=as.numeric(format(as.Date(mtbs@data$ALARM_D),"%d"))

fires=mtbs@data$GIS_ACR

frap = array(data = NA, dim = c(length(fires), 4)) #spi1,2,3,4,5,6,9,12,24,tasmax.tasmin,tas
frap[,1 ] = anni
frap[,2 ] = mesi
frap[,3 ] = giorni
frap[,4 ] = fires

writeMat(file.path(dir_fire, paste0("frap_baileys.mat")), frap = frap)
save(frap, file = file.path(dir_fire,"frap_baileys.RData"))

plot.ts(fires)
semilogy(anni,fires)
min(fires)


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



save(BA, file = file.path(dir_fire,"fire_1950_2020_monthly_baileys.RData"))
writeMat(file.path(dir_fire, paste("fire_1950_2020_monthly_baileys.mat",sep="")), BA=BA)

